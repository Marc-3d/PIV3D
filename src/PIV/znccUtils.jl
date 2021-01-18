# sumS and sumS2 are padded by 1 pixel at the start of each dimension, to avoid writting code to specifically deal with the first line and first column of the images. 

# 2D SEARCH AREA SUM-TABLES. 

function znccSumArrays!( search::AbstractArray{T,2}, sumS::A{U,2}, sumS2::A{U,2} 
                       ) where {T<:Real, U<:Real} 
    h, w = size(search);
    @inbounds for c in 2:w+1
        for r in 2:h+1
			 sumS[r,c] = search[h-r+2,w-c+2]   +  sumS[r-1,c] +  sumS[r,c-1] -  sumS[r-1,c-1]; 
            sumS2[r,c] = search[h-r+2,w-c+2]^2 + sumS2[r-1,c] + sumS2[r,c-1] - sumS2[r-1,c-1];
        end
    end
end

function znccSumValues( sumS::Array{T,2}, sumS2::Array{T,2},
                        row, col, ss, si ) where {T<:Real}
            
    pad = 1; 
    row = row + pad; 
    col = col + pad; 
    r0, r1 = max( 1, row - si[1] ), min( ss[1] + pad, row ); 
    c0, c1 = max( 1, col - si[2] ), min( ss[2] + pad, col ); 
        
    opS  = 0.0
    opS += sumS[ r1, c1 ]; 
    opS -= sumS[ r1, c0 ];
    opS -= sumS[ r0, c1 ];
    opS += sumS[ r0, c0 ];
        
    opS2  = 0.0
    opS2 += sumS2[ r1, c1 ]; 
    opS2 -= sumS2[ r1, c0 ];
    opS2 -= sumS2[ r0, c1 ];
    opS2 += sumS2[ r0, c0 ];
    
    return opS, opS2
end

# 3D SEARCH VOLUME SUM-TABLES.     

# optimized sumarea, reusing previous operations
function znccSumArrays!( SA::AbstractArray{U,3}, S::Array{T,3}, S2::Array{T,3} 
                       ) where {T<:Real,U<:Real} 

     h,  w,  d = size(SA);
    @inbounds for z in 2:d+1   
        for c in 2:w+1
            topSum  = 0.0; 
            topSum2 = 0.0; 
            for r in 2:h+1      
                SAVal  = SA[h-r+2,w-c+2,d-z+2]
                SAVal2 = SAVal*SAVal;
                    
                 S[r,c,z] = SAVal  +  S[r,c-1,z] +  S[r,c,z-1] -  S[r,c-1,z-1] +  topSum; 
                S2[r,c,z] = SAVal2 + S2[r,c-1,z] + S2[r,c,z-1] - S2[r,c-1,z-1] + topSum2;
                    
                 topSum  += SAVal; 
                topSum2  += SAVal2; 
            end
        end
    end
end


function znccSumValues( sumS::Array{T,3}, sumS2::Array{T,3}, 
                        row, col, zet, ss, si ) where {T<:Real}
    pad = 1; 
    row = row + pad; 
    col = col + pad; 
    zet = zet + pad; 
    r0, r1 = max( 1, row - si[1] ), min( ss[1] + pad, row ); 
    c0, c1 = max( 1, col - si[2] ), min( ss[2] + pad, col ); 
    z0, z1 = max( 1, zet - si[3] ), min( ss[3] + pad, zet );
    
    opS  = 0.0; 
    opS += sumS[ r1, c1, z1 ]
    opS -= sumS[ r0, c1, z1 ]
    opS -= sumS[ r1, c0, z1 ]
    opS -= sumS[ r1, c1, z0 ]
    opS += sumS[ r1, c0, z0 ]
    opS += sumS[ r0, c1, z0 ]
    opS += sumS[ r0, c0, z1 ]
    opS -= sumS[ r0, c0, z0 ]
        
    opS2  = 0.0; 
    opS2 += sumS2[ r1, c1, z1 ]
    opS2 -= sumS2[ r0, c1, z1 ]
    opS2 -= sumS2[ r1, c0, z1 ]
    opS2 -= sumS2[ r1, c1, z0 ]
    opS2 += sumS2[ r1, c0, z0 ]
    opS2 += sumS2[ r0, c1, z0 ]
    opS2 += sumS2[ r0, c0, z1 ]
    opS2 -= sumS2[ r0, c0, z0 ]
        
    return opS, opS2
end
