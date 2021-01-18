# Val{true} is passed at the last mp iteration and Val{false} during the previous iterations. 

function approxTranslation( cmatrix::Array{T,N}, alg::String, ::Val{true} ) where {T<:AbstractFloat,N}
    ccenter = div.( size( cmatrix ), 2 );
    peak1, maxVal = firstPeak( cmatrix );
    return subPixel( cmatrix, alg, peak1, maxVal ) .- ccenter;
end

function approxTranslation( cmatrix::Array{T,N}, alg::String, ::Val{false} ) where {T<:AbstractFloat,N}
    ccenter = div.( size( cmatrix ) , 2 );
    return firstPeak( cmatrix )[1] .- ccenter
end

function subPixel( cmatrix::Array{T,N}, alg::String, peak1::NTuple{N,I}, maxVal::T ) where {T<:AbstractFloat,N}
    
    refinement = zeros( T, N ); 

    if alg == "" || indexOnBorder( peak1, size( cmatrix ) )
        # no sub-pixel refinement
    else 
        neighs = getNeighbours( peak1, cmatrix );  
        
        # ensuring values are positive
        minVal = minimum( neighs );
        for idx in 1:length(neighs)
            neighs[idx] = 1 + neighs[idx] - minVal
        end
        
        if     alg == "gaussian"
            refinement = point3Gaussian( neighs, 1+maxVal-minVal );
        elseif alg == "centroid"
            refinement = centroid( neighs, maxVal );
        end 
         
        # nan occur when neighbours contains 0's, resulting in 0/0
        refinement = [ isnan(x) ? 0.0 : x for x in refinement ] 
    end

    return peak1 .+ refinement
end

# 3-POINT GAUSSIAN SUBPIXEL

# 3-point Gaussian subpixel 2D

function point3Gaussian( ns::Array{T,2}, mx::T ) where {T<:AbstractFloat}
    return point3_2D( log(ns[4]), log(ns[6]), log(ns[2]), log(ns[8]), log(mx) ); 
end

function point3_2D( up::T, down::T, left::T, right::T, mx::T ) where {T<:AbstractFloat}
    return [ ( up  - down )/( 2* up  - 4*mx + 2*down  ), 
             (left - right)/( 2*left - 4*mx + 2*right ) ];
end

# 3-point Gaussian subpixel 3D
function point3Gaussian( ns::Array{T,3}, mx::T ) where {T<:AbstractFloat} 
    return point3_3D( log(ns[13]), log(ns[15]), log(ns[11]), log(ns[17]), 
                      log(ns[ 5]), log(ns[23]), log(mx) ); 
end

function point3_3D( up::T, down::T, left::T, right::T, front::T, back::T, mx::T ) where {T<:AbstractFloat}
    return [ (  up  - down )/( 2*  up  - 4*mx + 2*down  ),
             ( left - right)/( 2* left - 4*mx + 2*right ),
             (front - back )/( 2*front - 4*mx + 2*back  ) ];
end

# CENTROID SUBPIXEL

# Centroid subpixel 2D 

function centroid( ns::Array{T,2}, mx::T ) where {T<:AbstractFloat}
	return centroid_2D( ns[4], ns[6], ns[2], ns[8], mx )
end

function centroid_2D( up::T, down::T, left::T, right::T, mx::T ) where {T<:AbstractFloat}
    return [ ( up  - down )/(  up  - mx + down  ), 
             (left - right)/( left - mx + right ) ];
end

# Centroid subpixel 3D

function centroid( ns::Array{T,3}, mx::T ) where {T<:AbstractFloat}
	return centroid_3D( ns[13], ns[15], ns[11], ns[17], ns[5], ns[23], mx )
end

function centroid_3D( up::T, down::T, left::T, right::T, front::T, back::T, mx::T ) where {T<:AbstractFloat}
    return [ (  up  - down )/(   up  - mx + down  ), 
             ( left - right)/(  left - mx + right ), 
			 (front - back )/( front - mx + back  ) ];
end
