""" PIV PARAMETER OBJECT """
mutable struct PIVParameters
         corr::Union{ZNCC,FFT}          # cross-correlation algorithm, "zncc" or "fft"
    interArea::Tuple{Int64,Int64,Int64} # interrogation area size(s) in pixels
 searchMargin::Tuple{Int64,Int64,Int64} # search margin(s) in pixels
      overlap::Tuple{Int64,Int64,Int64} # pixel distance(s) between contiguous interrogation areas
        mpass::Int64    # multi-pass depth
        width::Int64    # exclusion radius around 1st peak for finding 2nd peak
         peak::String   # one of "", "gaussian", "centroid"
     sigNoise::String   # one of "", "peak2peak"
      filtFun::Function # Filtering function to skip cross-correlating IA's without signal
    threshold::Float64  # skip cross correlation if filtFun < threshold
        mpFun::Function # scaling function for IA, SM and ovp during multipass
       mpReps::Array{Int64,1} # number of times each multipass iteration is repeated
        units::Tuple{DU,DU,DU,TU} # units for the 3 spatial dimensions, and temporal units
         args::Tuple{III,III,III,I,I,S,S,Function,F,Function,Array{I,1}}
end

# constructor
PIVParameters( corr::Union{ZNCC,FFT}, IA::UnI, SM::UnI, ovp::UnI, mpass::I, width::I, peak::S,
               sn::S, ffn::Function, th::F, mpfn::Function, mpReps::Array{I,1}, uns::UNTS 
             ) = PIVParameters( corr, toIII(IA), toIII(SM), toIII(ovp), mpass, width, peak, sn,
                                ffn, th, mpfn, mpReps, uns, ( toIII(IA), toIII(SM), toIII(ovp),
                                mpass, width, peak, sn, ffn, th, mpfn, mpReps ) 
                              )

function parseCorr( corr::String )
    lowercase(corr) == "zncc" &&  return ZNCC();
    lowercase(corr) == "fft"  &&  return  FFT();
end

# interrogation area size, search margin and overlap can be specified for each dimension,
# or a sigle value can be applied to all dimensions.
toIII( input::Int64 ) = ( input, input, input )
toIII( input::Tuple{Int64,Int64} ) = ( input[1], input[2], 0 )
toIII( input::Tuple{Int64,Int64, Int64} ) = input

# constructor with default values
function setPIVParameters(; interArea::UnI=32, searchMargin::UnI=0, overlap::UnI=0,
                            mpass::I=1, width::I=2, corr::S="fft",
                            peak::S="gaussian", sigNoise::S="",
                            filtFun=(x)->maxval(x), threshold::F=-1.0,
                            mpFun=(x,mp)->(x*mp), mpReps::Array{I,1}=[1,],
                            units=(1m,1m,1m,1second) )
    return PIVParameters( parseCorr(corr), interArea, searchMargin, overlap, mpass,
                          width, peak, sigNoise, filtFun, threshold, mpFun, mpReps, units )
end

function checkPIVParameters( dsize, N, params )
    !( N !== 2  ||  N !== 3 ) && error("Your input data isn't 2D nor 3D.")
    IA = params.interArea;
    if N == 2
         IA = params.interArea[1:2]
    end
    any( params.interArea .< 2 ) && error("interArea = 1 makes little sense.")
    any( IA .> dsize ) && error("interArea is bigger than the data")
    any( params.overlap .>= params.interArea ) && error("Overlap must be < IA size")

    if size(params.mpReps,1) !== params.mpass
        updateParameters!( params, :mpReps, ones( Int64, params.mpass ) )
    end
end

""" SYNTHETIC PARTICLES PARAMETER OBJECT """
mutable struct synthParameters
      n::Int64   # total number of particles, ignored when dens > 0
   dens::Int64   # particle density ( nºparticles/IA )
      d::Int64   # depth. If depth > 1, generates a volume
      w::Int64   # width of generated image/volume
      h::Int64   # height of generated image/volume
     i0::Float64 # max intensity of particles
     dt::Float64 # standard deviation of gaussian around each particles
     th::Float64 # standard deviation of gaussian around lightsheet
    err::Float64 # determines rendering radius around each pixel
      z::Int64   # light-sheet position. Default z = 0
  noise::Int64   # varible for controlling noise. Not really implemented.
    rad::Int64   # overwrites the rendering radius computed from err.
   mode::Bool    # switches between ( x, y, z ) or ( y, x, z ) coordinates
   args::Tuple{I,I,I,I,I,F,F,F,F,I,I,I,Bool}
end

# constructor
synthParameters( n::Int64, dens::Int64, d::Int64, w::Int64, h::Int64, i0::Float64, dt::Float64,
                 th::Float64, err::Float64, z::Int64, noise::Int64, rad::Int64, mode::Bool 
               ) = synthParameters( n, dens, d, w, h, i0, dt, th, err, z, noise, rad, mode, 
                                   (n, dens, d, w, h, i0, dt, th, err, z, noise, rad, mode)
                                  )

# constructor with default values
function setSynthParameters(; n=0::I, dens=0::I, d=1::I, w=300::I, h=300::I,
                              z=0::I, i0=255.0::F, dt=3.0::F, err=0.1::F,
                              th=10.0::F, noise=0::I, rad=0::I, mode="whd"::S )

    return synthParameters( n, dens, d, w, h, i0, dt, th, err, z, noise, rad, mode=="whd" )
end

""" TRANSFORMATION PARAMETERS """
mutable struct transformParameters{N,M}
      kind::Symbol            # currently ony :translation is supported
     means::NTuple{N,Float64} # mean value of transform in x, y [ and z ] if N = 3
      vars::NTuple{N,Float64} # variance component in covariance matrix
 cov_ratio::NTuple{M,Float64} # covariance components in covariance matrix. -1 < 1
      args::Tuple{Symbol,NTuple{N,Float64},NTuple{N,Float64},NTuple{M,Float64}}
end

# constructor
transformParameters( kind::Symbol, means::NTuple{N,Float64}, vars::NTuple{N,Float64},
                     covs::NTuple{M,Float64}  
                   ) where {N,M} = transformParameters{N,M}( kind, means, vars, covs, 
                                                            (kind, means, vars, covs) 
                                                           )

# constructor with default values
function setTransformParameters(; kind=:translation, means=(1.0, 1.0), vars=(eps(Float64),
								  eps(Float64)), cov_ratio=(0.0,) )

    return transformParameters( kind, means, vars, cov_ratio );
end

""" EVALUATION PARAMETER OBJECT """
mutable struct metaParameters
 metaloop::Int64    # determines the
 variable::Symbol   # what variable will be changed each iteration.
      min::Float64  # minimum value of changing variable
      max::Float64  # maximum value of changing variable
  repeats::Int64    # number of measurements for computing bias and error
     args::Tuple{Int64,Symbol,Float64,Float64,Int64}
end

metaParameters( metaloop::Int64, variable::Symbol, min::Float64, max::Float64, repeats::Int64
              ) = metaParameters( metaloop, variable, min, max, repeats, 
                                 (metaloop, variable, min, max, repeats)
                                )

# short type aliases, to make function signatures shorter
PP   = PIVParameters;
SP   = synthParameters;
TP2D = transformParameters{2,1}
TP3D = transformParameters{3,3}
TP   = Union{TP2D,TP3D}
MP   = metaParameters;



"""
    Creates a set of parameters for directing PIV evaluation. Accepted parameters are:
        variable (  Symbol ): a symbol indicating what variable should be evaluated. Can be any parameters
                              in synthParameters, PIVParameters or transformParameters.
        metaloop (  Int64  ): number of iterations
            min  ( Float64 ): min value of changing variable
            max  ( Float64 ): max value of changing variable
        repeats  (  Int64  ): number of repetitions
"""
function setMetaParameters( ; metaloop=0::I, var=Symbol("")::Symbol, min=0.0::F,  max=0.0::F, repeats=1::I )

    return metaParameters( metaloop, var, min,  max, repeats )
end

# Auxiliary functions

to2DTransform( t2D::TP2D ) = t2D

function to2DTransform( t3D::TP3D )
    means     = (  t3D.means[1], t3D.means[2] );
    vars      = (  t3D.vars[1] ,  t3D.vars[2] );
    cov_ratio = (  t3D.cov_ratio[1], );

    return transformParameters{2,1}( t3D.kind, means, vars, cov_ratio );
end

function to3DTransform( t2D::TP2D )
    means     = (   t2D.means[1]  ,   t2D.means[2]  , (t2D.means[1] + t2D.means[2])/2 );
    vars      = (    t2D.vars[1]  ,    t2D.vars[2]  , ( t2D.vars[1] +  t2D.vars[2])/2 );
    cov_ratio = ( t2D.cov_ratio[1], t2D.cov_ratio[1], t2D.cov_ratio[1] );

    return transformParameters{3,3}( t2D.kind, means, vars, cov_ratio );
end

to3DTransform( t3D::TP3D ) = t3D


# function to update parameters, and also update the .args tuple, which was created to pass parameters
# by tuple deconstruction "parameter.args..."

function updateParameters!( params, field, newVal )
    hasfield = false;
    fieldidx =     0;
    for e in fieldnames( typeof(params) )
        fieldidx += 1;
        if ( e == field )
            hasfield = true;
            setproperty!( params, e, newVal )
            break
        end
    end

    ( !hasfield ) && return nothing;

    if typeof(params) == PIVParameters
        fieldidx -= 1 # ignore corr parameter
    end
    newargs = [ arg for arg in params.args ]
    newargs[ fieldidx ] = newVal;
    setproperty!( params, :args, Tuple(newargs) )
end
