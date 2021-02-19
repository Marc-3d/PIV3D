module PIV3D

using FFTW, FileIO, LIBTIFF, Statistics

export PIV

# Shortened type names to keep function signatures short ( when possible one-liners )
include("units.jl")

F      = Float64;
I      = Int64;
A      = Array;
S      = String;
C      = Complex;
II     = Tuple{Int64,Int64};
III    = Tuple{Int64,Int64,Int64};
UnI    = Union{I,II,III}
SySySy = Tuple{Symbol,Symbol,Symbol}
DU     = DistanceUnits
TU     = TimeUnits
UNTS   = Tuple{DU,DU,DU,TU}

struct ZNCC end
struct FFT  end

include("parameters.jl")
include("pivIO.jl")
include("utils.jl")

include("PIV/gaussianSubpixel.jl")
include("PIV/signalToNoise.jl")
include("PIV/ccrZNCC.jl")
include("PIV/ccrFFT.jl")

include("postprocessing.jl")


function PIV( data1::Array{T,N}, data2::Array{T,N}, params::PIVParameters; checkParams=true, units=nothing ) where {T<:Real,N}

    checkParams && checkPIVParameters( size(data1), N, params )

    #Continuing to PIV analysis
    _PIV( data1, data2, params )
end

""" Calling 2D PIV analysis """
function _PIV( img1::A{T,2}, img2::A{T,2}, p::PIVParameters; units=nothing ) where {T<:Real}

    u, v, sn = PIV_2D( p.corr, img1, img2, p.args... );

    if units !== nothing
        ratios = changeUnits.(  p.units[1:2], units[1:2] )./changeUnits( p.units, units[4]  )
        for e in 1:length(u)
            u[e] = u[e]*ratios[1]
            v[e] = v[e]*ratios[2]
        end
    end
    return u, v, sn;
end

""" Calling 3D PIV analysis """
function _PIV( vol1::A{T,3}, vol2::A{T,3}, p::PIVParameters; units=nothing ) where {T<:Real}

    u, v, w, sn = PIV_3D( p.corr, vol1, vol2, p.args... );

    if units !== nothing
        ratios = changeUnits.( p.units, units[1:3] )./changeUnits( p.units, units[4] )
        for e in 1:length(u)
            u[e] = u[e]*ratios[1]
            v[e] = v[e]*ratios[2]
            w[e] = w[e]*ratios[3]
        end
    end

    return u, v, w, sn;
end

"""
    Automated PIV analysis between each pair of consecutive timepoints from a time-sequence data.

    This function assumes that the files follow a consistent naming scheme, where the name of each
    file contains a varying index indicating its time index. Aside from the varying index, the rest
    of the file name should be the same between all files. For example:

        result_001_monday, result_002_monday, result_003_monday, etc..

    To run a sequential PIV analysis on the previous example the following arguments should be give:

        sequencePIV( pp, "result_", "_monday", 1, 3, digits=3 )

    where pp is a PIV parameter object, "result_" and "_monday" are the invariant parts in each file name
    before and after the varying index, 1 is the first index, 3 is the last index to analyze, and digits=3
    indicates the number of digits of the varying index. 0's are prepended to reach "digits".

"""
function sequencePIV( PIVparams::PIVParameters,
                      before_id::String, after_id::String,
                      id0::Integer, idStop::Integer;
                      digits::Integer=1,
                      path::String=pwd(),
                      subsample::Int64=1,
                      typ=Float32 )

    n = idStop - id0

    U = Array{ Array{Float32,3}, 1 }( undef, n );
    V = Array{ Array{Float32,3}, 1 }( undef, n );
    W = Array{ Array{Float32,3}, 1 }( undef, n );

    if path[end] !== '/'
        path *= '/'
    end

    for idx in id0:idStop-1

        idx1  = numberDigits(  idx , maxdigits=digits );
        idx2  = numberDigits( idx+1, maxdigits=digits );
        file1 = before_id * idx1 * after_id;
        file2 = before_id * idx2 * after_id;
        vol1  = PIVload( path*file1, sample=subsample, typ=typ );
        vol2  = PIVload( path*file2, sample=subsample, typ=typ );

        println( "Performing PIV between:" )
        println( "\t$(path*file1)" )
        println( "\t$(path*file2)" )

        u, v, w = PIV( vol1, vol2, PIVparams );

        U[1+idx-id0] = u;
        V[1+idx-id0] = v;
        W[1+idx-id0] = w;
    end

    return U, V, W
end


# projectionStackPIVOFF in smbnotes
"""
    PIV3D.projectionStackPIV( PIVParams, filename )

    Computes 2D PIV between stack[idx] and stack[idx+1], for all consecutive
    projections in the stack.
"""
function projectionStackPIV( PIVParams, fn::String; threshold=-1, lo=1, hi=-1, reps=Int64[], typ=nothing )

    stack = PIVload( fn, typ=typ );
    size( stack, 3 ) < 2 && error("Stack only contains one image, no PIV possible")

    lo = ( lo < size( stack, 3 ) ) ? lo :  1;
    hi = ( hi == -1 ) ? size(stack,3)-1 : hi;

    i1 = stack[ :, :,  lo  ];
    i2 = stack[ :, :, lo+1 ];

    u, v, _  = PIV( i1, i2, PIVParams );

    U = zeros( size(u,1), size(u,2), hi - lo + 2 );
    V = zeros( size(v,1), size(v,2), hi - lo + 2 );

    U[:,:,1] .= u;
    V[:,:,1] .= v;

    for idx in lo+1:hi
        println("Analyzing $(idx-lo)/$(hi-lo)");

        i1 = stack[ :, :,  idx  ];
        i2 = stack[ :, :, idx+1 ];

        u, v, _ = PIV( i1, i2 , PIVParams );

        U[:,:,idx-lo+1] .= u;
        V[:,:,idx-lo+1] .= v;
    end

    return U, V
end

function projectionStackPIV( PIVParams, stack::Array{T,3};
                             lo=1, hi=-1, reps=Int64[] ) where {T<:Real}

    lo = ( lo < size( stack, 3 ) ) ? lo :  1;
    hi = ( hi == -1 ) ? size(stack,3)-1 : hi;

    i1 = stack[ :, :,  lo  ];
    i2 = stack[ :, :, lo+1 ];

    u, v, _  = PIV( i1, i2, PIVParams );

    U = zeros( size(u,1), size(u,2), hi - lo + 2 );
    V = zeros( size(v,1), size(v,2), hi - lo + 2 );

    U[:,:,1] .= u;
    V[:,:,1] .= v;

    for idx in lo+1:hi
        println("Analyzing $(idx-lo)/$(hi-lo)");

        i1 = stack[ :, :,  idx  ];
        i2 = stack[ :, :, idx+1 ];

        u, v, _ = PIV( i1, i2 , PIVParams );

        U[:,:,idx-lo+1] .= u;
        V[:,:,idx-lo+1] .= v;
    end

    return U, V
end

# Evaluation functions and plotting to be used in Jupyter notebooks

using DataFrames, Gadfly, Compose, Distributions

include("EvaluationOnSyntheticData/syntheticParticles.jl")
include("EvaluationOnSyntheticData/transformations.jl")
include("EvaluationOnSyntheticData/evaluation.jl")

end
