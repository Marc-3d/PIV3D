export generateParticles, renderSlice, normalize, renderVolume, volumeToGif

"""
    computeIntensityGaussian( params::synthParameters, x::Int64, y::Int64, z::Int64, x0::Real, y0::Real, z0::Real )

    Computes the intensity contribution of particles ( x0, y0, z0 ) to pixel ( x, y ) using the following formula
    (M.Raffael 1998 Particle Image Velocimetry):

            I(x,y) = MaxI(z0)*exp( -1 * ( (x - x0)² + (y - y0)² )/(dt²/8) )

            where:
                     MaxI(z0) = i0*exp( -1 * (z - z0)²/(th²/8) )

    If the lightsheet is infinitely thick ( th = Inf64 ) or (z0 = z) then MaxI(z0) = i0.
"""

function computeIntensityGaussian( params::SP, x::I, y::I, z::I, x0::T, y0::T, z0::T ) where {T<:Real}
    return computeIntensityGaussian( params.i0, params.dt, params.th, x, y, z , x0, y0, z0 )
end

function computeIntensityGaussian( i0::F, dt::F, th::F, x::I, y::I, z::I, x0::T, y0::T, z0::T ) where {T<:Real}

    return i0*exp( -8*( (z - z0)^2/th^2 + ((x - x0)^2 + (y - y0)^2)/dt^2 ) )
end

#=
    Instead of calculating the Gaussian intensity of each particle over the whole image, it is more efficient to
    limit the area over which we compute the gaussian for each particle. This method calculates the radius ( in
    pixels, aka Int64 ) at which the Gaussian reaches a certain intensity, named "error". The formula is obtained
    by considering a 1D gaussian, and the distance "d" at which the gaussian equals "error".

    I( d ) = i0 * e^( ( - d² ) / ( 1/8 * dt²) ) = error
           =>          ln(i0) - d² / ( dt²/8 )  = ln(error)
           =>                                 d = sqrt( dt²( ln(i0) - ln(error) )/8 )

    Same goes for the z-axis. I think the zRadius needs to include the effect of dt, not sure.
=#

calculateGaussianRadius( params::SP ) = calculateGaussianRadius( params.err, params.i0, params.dt, params.rad )

function calculateGaussianRadius( error::F, i0::F, dt::F, rad::I )
    return ( rad > 0 ) ? rad : round( Int64, sqrt( dt*dt*( log(i0) - log(error) )/8 ) )
end

calculateGaussianZRadius( params::SP ) = calculateGaussianZRadius( params.err, params.i0, params.th, params.d )

function calculateGaussianZRadius( error::F, i0::F, th::F, depth::I )
    return ( th == Inf ) ? depth : round( Int64, sqrt( th*th*( log(i0) - log(error) )/8 ) )
end

#=
    syntheticParticles.generateParticles( synParams, [pivparams=pivParams] );
    ---

    Generates the coordinates of the synthetic particles, with the data type ::Array{ Tuple{ Float64,Float64,Float64 }, 1 }.
    Particle coordinates are Float's because PIVLab does the same, and it increases overall accuracy assessment at decimal
    translations.

    if ( synParams.dens > 0 )

        generateDensityParticles( args ) is called. synParams.dens determines the number of particles in each
        interrogation area. The size of the interrogation area is extracted from a pivParameters object, which is
        passed as an optional argument.

        example:
                    synParams = syntheticParticles.setSynthParameters();
                    pivParams = nopad_PIV.setPIVParameters()
                    coords    = syntheticParticles.generateParticles( synParams, pivparams=pivParams )

    else

        generateParticles( args ) is called. synParams.n is used to determine the total number of particles in the
        generated volume.

        example:
                    synParams = syntheticParticles.setSynthParameters();
                    coords    = syntheticParticles.generateParticles( synParams ).


    Returned particles are sorted in the Z-axis.
=#

newParticle( ::Val{true} , y0, x0, z0, ys, xs, zs ) = ( x0 + rand()*xs, y0 + rand()*ys, z0 + rand()*zs );
newParticle( ::Val{false}, y0, x0, z0, ys, xs, zs ) = ( y0 + rand()*ys, x0 + rand()*xs, z0 + rand()*zs );

function generateParticles( p::SP; pivparams=setPIVParameters()::PIVParameters )
    if p.dens > 0
        return generateDensityParticles( p.dens, p.w, p.h, p.d, pivparams.interArea, p.mode );
    else
        return generateTotalNumParticles( p.n, p.w, p.h, p.d, p.mode )
    end
end

function generateTotalNumParticles( n::I, w::I, h::I, d::I, mode::Bool )

    particles = [ newParticle( Val(mode), 0, 0, 0, h, w, d ) for x in 1:n ];
    return sort!( particles, by = x -> x[3] )
end

function generateDensityParticles( dens::I, w::I, h::I, d::I, IA::III, mode::Bool )

    particles  = Array{ Tuple{Float64,Float64,Float64}, 1}(undef,0);
    Y0, X0, Z0 = min.( IA, ( h, w, d ) );
    thicknessZ = ( d == 1 ) ? 1 : IA[3];


    for z_max in Z0:IA[3]:d
        for x_max in X0:IA[2]:w
            for y_max in Y0:IA[1]:h
                for p in 1:dens
                    push!( particles, newParticle( Val(mode), y_max-Y0, x_max-X0, z_max-Z0, IA[1], IA[2], thicknessZ ) );
                end
            end
        end
    end

    return sort!( particles, by = x -> x[3] );
end


function renderParticles( particles::Array{Tuple{T,T,T},1}, params::SP ) where {T<:Real}
    if params.d == 1
        return renderSlice( particles, params );
    else
        return renderVolume( particles, params );
    end
end

renderSlice( particles::A{Tuple{T,T,T},1}, p::SP ) where {T<:Real} = renderSlice( particles, p.w, p.h, p.d, p.i0, p.dt, p.th, p.z, p.err, p.noise, p.rad, p.mode )


function renderSlice( particles::A{Tuple{T,T,T},1}, w::I, h::I, d::I, i0::F, dt::F, th::F, z::I,
                      err::F, noise::I, rad::I, mode::Bool ) where {T<:Real}

    img = zeros( Float64, ( h, w ) );
    rz  = calculateGaussianZRadius( err, i0, th,   d );
    r   =  calculateGaussianRadius( err, i0, dt, rad );

    for p in particles
        if abs( p[3] - z ) <= rz

            xmin, xmax = max( 1, floor( Int, p[1] - r ) ), min( w, ceil( Int, p[1] + r ) );
            ymin, ymax = max( 1, floor( Int, p[2] - r ) ), min( h, ceil( Int, p[2] + r ) );

            for x in xmin:xmax
                for y in ymin:ymax
                    @inbounds img[ y, x ] += computeIntensityGaussian( i0, dt, th, x, y, z, p[1], p[2], p[3] )
                end
            end
        end
    end

    # Allow for noise. Only Poisson noise is implemented by now.
    if ( noise > 0 )
        poi = Poisson( noise );
        for e in eachindex( img )
            img[e] += rand( poi, 1 )[1]
        end
    end

    return img
end


renderVolume( particles::A{Tuple{T,T,T},1}, p::SP ) where {T<:Real} = renderVolume( particles, p.w, p.h, p.d, p.i0, p.dt, p.th, p.err, p.rad );

function renderVolume( particles::A{Tuple{T,T,T},1}, w::I, h::I, d::I, i0::F, dt::F, th::F, err::F, rad::I ) where { T<:Real }

    out = zeros( Float64, h, w, d );
    r   = calculateGaussianRadius( err, i0, dt, rad );

    #println( "[h,w,d] = ", [h,w,d], " radius = ", r );

    for p in particles

        p_lo, p_hi = floor.( Int, p .- r ), floor.( Int, p .+ r );
        ymin, ymax = max( 1, p_lo[1] ), min( h, p_hi[1] );
        xmin, xmax = max( 1, p_lo[2] ), min( w, p_hi[2] );
        zmin, zmax = max( 1, p_lo[3] ), min( d, p_hi[3] );

        for z in zmin:zmax
            for x in xmin:xmax
                for y in ymin:ymax
                    @inbounds out[ y, x, z ] += computeIntensityGaussian( i0, dt, dt, y, x, z, p[1], p[2], p[3] )
                end
            end
        end
    end

    return out
end


"""
    normalize( img::Array{Float64,2}, range=(0,0)::Tuple{R, R} )
    ---
    Normalizes the intensities values between 0.0 and 255.0.

    Normalization is not automatically applied because it can change the intensity values of pixels when there is
    out-of-frame loss of the particles with the maximum intensities. Therefore, when normalizing the original and
    tranformed image, the minimum and maximum values of the original image should be used in the normalization
    of the transformed image. Example:

                     slice     = syntheticParticles.renderSlice( coords, synParams );
                     slice_t   = syntheticParticles.renderSlice( transformed_coords, synParams );

                     min, max  = extrema(slice)
                     slice      = syntheticParticles.normalize( slice );
                     slice_t    = syntheticParticles.normalize( slice_t , ( min, max ) );
"""
function normalize( img::AbstractArray{T,N} ) where {T<:Real,N}

    min, max = extrema( img );
    return normalize( img, ( min, max ) )
end

function normalize( img::AbstractArray{T,N}, range::Tuple{T,T} ) where {T<:Real,N}

    min, max = range;
    norm_img = zeros( eltype(img), size(img) );
    diff     = max - min;
    @simd for i in 1:length(norm_img)
       @inbounds norm_img[i] = (img[i] - min)/diff
    end
    return norm_img
end
