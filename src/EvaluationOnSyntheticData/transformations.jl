export applyTransformation


# Sample transformation, right now only translations are possible

function transformation( params::transformParameters ) 
    params.kind == :translation && return translation( params ); 
                                   return false; 
end

"""
    transformations.translation( transformParams::transformParameters )
    ---
    
    Samples a translation from a multivariate normal distribution. Using this approach, 
    we can test our algorithm on translation of different strengths and directions, and 
    test the accuracy of our algorithm on more heterogeneous conditions.

    NOTE: the variance-covariance matrix passed to the multivariate normal distribution 
    needs to be hermitian and positive definite. This means that, when the users want, 0 
    variance in their samples, they should pass: vars=[ eps(Float64), eps(Float64) ], which 
    is the smallest non-zero value a Float64 can take. This reflects that the product of 
    variances is always greater than the product of coviarances. 
"""
         translation( params::transformParameters ) = translation( params.means, params.vars, params.cov_ratio )

function translation( means::NTuple{N,Float64}, vars::NTuple{N,Float64}, cov_ratio::NTuple{M,Float64} ) where {N,M} 
        
    dist = Distributions.MvNormal( [ x for x in means ], covarmatrix( vars, cov_ratio ) ); 
    return Tuple( rand( dist, 1 ) ); 
end

function covarmatrix( vars::NTuple{2,Float64}, cov_ratio::Tuple{Float64} ) # 2D covariance matrix
    
    cov = sqrt( vars[1]*vars[2] )*cov_ratio[1]; 
    return [ vars[1] cov; cov vars[2] ]
end

function covarmatrix( vars::NTuple{3,Float64}, cov_ratio::NTuple{3,Float64} )  # 3D covariance matrix
    
    cov1 = sqrt( vars[1]*vars[2] )*cov_ratio[1]; 
    cov2 = sqrt( vars[2]*vars[3] )*cov_ratio[2]; 
    cov3 = sqrt( vars[1]*vars[3] )*cov_ratio[3]; 
    return [ vars[1] cov1 cov3; cov1 vars[2] cov2; cov3 cov2 vars[3] ]
end

# Return 3D particles after applying a transformation

function applyTransformation( particles::Array{NTuple{3,T},1}, tp::transformParameters ) where {T<:Real}
    t = transformation( tp )
    return applyTransformation( particles, t )
end

# Apply 2D transform to 3D particles
function applyTransformation( particles::Array{NTuple{3,T},1}, t::NTuple{2,Float64} ) where {T<:Real}
    return [ ( x[1] + t[1], x[2] + t[2], x[3] + 0.0 ) for x in particles ]
end

# Apply 3D transform to 3D particles 
function applyTransformation( particles::Array{NTuple{3,T},1}, t::NTuple{3,Float64} ) where {T<:Real}
    return [ ( x[1] + t[1], x[2] + t[2], x[3] + t[3] ) for x in particles ]
end

