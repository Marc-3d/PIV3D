"""
   SIGNAL TO NOISE MEASUREMENT
"""
function SNRatio( corr::Array{T,N}, method::S, width::I ) where {T<:Real,N}
        
    max2       = 0.0; 
    h, w, d    = size(corr)
    idx1, max1 = firstPeak( corr )
    
    if     method == "peak2peak"
        idx2, max2 = secondPeak( corr, width, idx1 )
    elseif method == "peak2mean"
        max2 = mean(corr)
    else
        throw(ErrorException("Error: use 'peak2peak' or 'peak2mean'.\n"))
    end
    
    # if max2 == 0, sig2noise = inf. Check this at return
    sig2noise = max2 != 0.0 ? max1 / max2 : 0.0
    return sig2Noise; 
end
