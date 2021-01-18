#=
	2D IMPLEMENTATION
=#

# crosscorr( f, g ) = iFFT( conj( FFT( f ) ) .* FFT( g ) )
function crossCorrelation!( ::FFT, f::A{C{T},2}, g::A{C{T},2}, corr::A{T,2}, plan, iplan 
                          ) where {T<:AbstractFloat} 
    plan * f;
	plan * g;
    @inbounds @simd for e in 1:length(f)
        f[e] = conj( f[e] ) * g[e]; 
    end
    
	iplan * f; 
	@inbounds @simd for e in 1:length(f) 
		corr[e] = real( f[e] ) 
	end
end

function PIV_2D( ::FFT, img1::A{T,2}, img2::A{T,2},
                        IA::III, SM::III, overlap::III, mpass::I, width::I,
                        peak::S, sigNoise::S, filtFun::Function, threshold::F, 
                        mpFun::Function, reps::Array{I,1}; 
                        corrType=Float32 
               ) where {T<:Real} 
        
	# Calculating vector field size.
    h, w   = size( img1 ); 
    step   = IA  .-  overlap; 
    VFsize = 1 .+ floor.( Int64, ((h,w,0).-IA)./step ); 
    U  = zeros( Float64, VFsize[1:2] ); 
    V  = zeros( Float64, VFsize[1:2] ); 
    SN = zeros( Float64, VFsize[1:2] );

	ignoreSN = sigNoise == ""; 

    @inbounds for mp in mpass:-1:1
                  
        # Scaling IA, SM and step (overlap) to multi-pass.
        IA_mp   = mpFun.(  IA , mp ); 
        SM_mp   = mpFun.(  SM , mp ); 
        step_mp = mpFun.( step, mp ); 
        last_mp = ( mp == 1 ); 
        
        ssize   = IA_mp .+ 2 .* SM_mp
        search  = zeros( T, ssize[1:2] ); 
                        
        # Initializing FFT crossCorrelation variables. 
        csize   = 2 .* ( IA_mp .+ SM_mp )
        cmatrix = zeros( corrType, csize[1:2] );
        shifted = zeros( corrType, csize[1:2] ); 
        pads    = zeros( Complex{corrType}, csize[1:2] );
        padi    = zeros( Complex{corrType}, csize[1:2] ); 
        plan    = FFTW.plan_fft!(  pads ); 
        iplan   = FFTW.plan_ifft!( pads ); 
        shifts  = div.( csize, 2 )[1:2] .- SM_mp[1:2] .- 1; 
            
        # Constructing grid coordinates for each interrogation area.
        grid = constructGrid( (h, w, IA_mp[3]), IA_mp, step_mp, mp ) 
            
        for g in grid
            #  g[1],  g[2] are start and stop rows of interrogation area
            #  g[3],  g[4] are start and stop columns of interrogation area
            #  g[7],  g[8] are start and stop rows inside UV arrays
            #  g[9], g[10] are start and stop columns inside UV arrays
                
            interr = view( img1, g[1]:g[2], g[3]:g[4] )

            if threshold > 0 && filtFun( interr ) < threshold 
                continue;
            end
                
            # 1-. Shifting the searchArea ranges by previously calculated displacements. 
            offH = round(Int64, U[ g[7], g[9] ]);
            offW = round(Int64, V[ g[7], g[9] ]);
                
            # 1.2-. Search volume coordinates after shifting.
            scoords, so = searchGrid( g, SM, offH, offW, 0, h, w, 0 )
                
            setTo0!( search )
            putWithinSearch!( search, img2, so, scoords )
                
            setTo0!( pads )
            setTo0!( padi )
            putWithinPadded!( padi, interr, 0.0 ); 
            putWithinPadded!( pads, search, 0.0 );
                        
            # 2-. cross-correlation
            crossCorrelation!( FFT(), pads, padi, cmatrix, plan, iplan )
            Base.circshift!( shifted, cmatrix, shifts )
                
            # 3-. calculating displacement
            ( r, c, ) = approxTranslation( shifted, peak, Val(last_mp) )
			if ( !ignoreSN && last_mp )
            	SN[ g[7], g[9] ] = SNRatio( cmatrix, sigNoise, width )
			end

            # 4-. Updating U, V matrices 
            t1 =  g[7]:g[8]; 
            t2 =  g[9]:g[10]; 
            U[ t1, t2 ] .= offH - r;
            V[ t1, t2 ] .= offW - c;
            
        end
    end # Multi-pass loop

    return U, V, SN
end

#= 
    3D IMPLEMENTATION
=#

function crossCorrelation!( ::FFT, f::A{C{T},3}, g::A{C{T},3}, corr::A{T,3}, plan, iplan 
                          ) where {T<:AbstractFloat} 

    plan * f;
    plan * g;
    @inbounds @simd for e in 1:length(f)
        f[e] = conj( f[e] ) * g[e]; 
    end

    iplan * f; 
    @inbounds @simd for e in 1:length(f) 
        corr[e] = real( f[e] ) 
    end
end

function PIV_3D( ::FFT, vol1::A{T,3}, vol2::A{T,3},
                        IA::III, SM::III, overlap::III, mpass::I, width::I,
                        peak::S, sigNoise::S, filtFun::Function, threshold::F, 
                        mpFun::Function, reps::Array{I,1}; corrType=Float32
               ) where {T<:Real} 
    

	# Calculating vector field size.
    h, w, d = size( vol1 );
    step   = IA .- overlap; 
    VFsize = 1 .+ floor.( Int64, (size(vol1) .- IA)./step ); 
    U  = zeros( corrType, VFsize ); 
    V  = zeros( corrType, VFsize ); 
    W  = zeros( corrType, VFsize );
    SN = zeros( corrType, VFsize );

	ignoreSN = sigNoise == ""; 

    @inbounds for mp in mpass:-1:1
                
        # Scaling IA, SM and step (overlap) to multi-pass iteration.
        IA_mp   = mpFun.(  IA , mp ); 
        SM_mp   = mpFun.(  SM , mp ); 
        step_mp = mpFun.( step, mp ); 
        last_mp = ( mp == 1 ); 
            
        # Initialize cross-correlation variables once per multi-pass iteration
        csize   = 2 .* ( IA_mp .+ SM_mp )
        shifts  = div.( csize, 2 ) .- SM_mp .- 1; 
        cmatrix = zeros( corrType, csize );
        shifted = zeros( corrType, csize ); 
        padi    = zeros( Complex{corrType}, csize ); 
        pads    = zeros( Complex{corrType}, csize );
        plan    = FFTW.plan_fft!(  padi ); 
        iplan   = FFTW.plan_bfft!( padi );
                
        # Constructing grid coordinates for each interrogation volume. 
        grid = constructGrid( size(vol1), IA_mp, step_mp, mp ) 
            
        for g in grid
            #  g[1],  g[2] are start and stop rows of interrogation volume
            #  g[3],  g[4] are start and stop columns of interrogation volume
            #  g[5],  g[6] are start and stop depths of interrogation volume
            #  g[7],  g[8] are start and stop rows inside UVW arrays
            #  g[9], g[10] are start and stop columns inside UVW arrays
            # g[11], g[12] are start and stop depths insie UVW arrays
            
            interr = view( vol1, g[1]:g[2], g[3]:g[4], g[5]:g[6] ) 

            if threshold > 0 && filtFun(interr) < threshold
                continue;
            end
                
            # 1-. Shifting the searchArea ranges by previously calculated displacements. 
            offH = round(Int64, U[ g[7], g[9], g[11] ]);
            offW = round(Int64, V[ g[7], g[9], g[11] ]);
            offD = round(Int64, W[ g[7], g[9], g[11] ]);
                
            # 1.2-. Search volume coordinates after shifting.
            s, so = searchGrid( g, SM, offH, offW, offD, h, w, d )
   
            setTo0!( padi, pads ); 
            putWithinPadded!( padi, interr, 0.0 )
            putWithinPadded!( pads, vol2  , 0.0, so, (s[1],s[2],s[3],s[4],s[5],s[6]) );
             
            # 2-. Computing FFT cross correlation
            crossCorrelation!( FFT(), pads, padi, cmatrix, plan, iplan );
            Base.circshift!( shifted, cmatrix, shifts );

            # 3-. Calculating displacements
            ( r, c, z ) = approxTranslation( shifted, peak, Val(last_mp) )
 			if ( !ignoreSN && last_mp )
            	SN[ g[7], g[9], g[11] ] = SNRatio( cmatrix, sigNoise, width )
			end

            # 4-. Updating U, V, W matrices 
            t1 = g[7]:g[8]; 
            t2 = g[9]:g[10]; 
            t3 = g[11]:g[12]; 
            U[ t1, t2, t3 ] .= offH - r;
            V[ t1, t2, t3 ] .= offW - c;
            W[ t1, t2, t3 ] .= offD - z; 
        end                   
    end # Multi-pass loop
    
    return U, V, W, SN
end



""" FFT out of place CROSS CORRELATION """
function crossCorrelation(::FFT, f::A{T,N}, g::A{T,N}; shift=false ) where {T<:Real,N} 
    
    corrSize = size( f ) .+ size( g ); 

    corr = zeros( Float64, corrSize ); 
    padf = zeros( Complex{Float64}, corrSize ); 
    padg = zeros( Complex{Float64}, corrSize ); 

    padf[ 1:size(f,1), 1:size(f,2), 1:size(f,3) ] .= complex.(f, 0.0 ); 
    padg[ 1:size(g,1), 1:size(g,2), 1:size(g,3) ] .= complex.(g, 0.0 ); 

    plan  =  plan_fft!( padf );
    iplan = plan_ifft!( padf ); 

    plan * padf;
    plan * padg;
    @inbounds @simd for e in 1:length(padf)
        padf[e] = conj( padf[e] ) * padg[e]; 
    end

    iplan * padf; 
    @inbounds @simd for e in 1:length(padf) 
        corr[e] = real( padf[e] ) 
    end
   
    if shift
        shifted = copy( corr ); 
        shifts = div.( corrSize, 2 ) .- div.( size(f) .- size(g), 2 ) .- 1; 
        Base.circshift!( corr, shifted, shifts );
    end
    
   return corr
end
