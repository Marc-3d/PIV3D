

# Used to load sequential data where the time index is preceeded by 0's: 0004
function numberDigits( number; maxdigits=4 )
   ndigs = 0
   for d in 0:maxdigits
       ndigs += number ÷ 10^d > 0
   end
   nzeros = ( maxdigits - ndigs );
   return ( nzeros < 0 ) ? string(number) : "0"^(nzeros)*string(number)
end



function vectorFieldSize( datasize, pp::PIVParameters )
    step = pp.interArea .- pp.overlap
    return 1 .+ floor.( Int64, ( datasize .- pp.interArea )./step )
end



function constructGrid( datasize::III, IA::III, step::III, mp::I )

    dims  = 1 .+ floor.( Int64, ( datasize .- IA )./step )
    grid  = Array{NTuple{12,Int64}}( undef, prod( dims ) )
    count = 1;

    # Each square/cube subdivision of the input images/volumes includes the image/volume
    # coordinates of the subdivision, and the coordinates of computed displacement
    # vector in the final vector field ( U, V, W ).
    rdelta = [ step[1], step[1], mp, mp ]
    cdelta = [ step[2], step[2], mp, mp ]
    ddelta = [ step[3], step[3], mp, mp ]

    @inbounds begin
    # ds = [ minDepthData, maxDepthData, minDepthUVW, maxDepthUVW ]
    ds = [ 1, IA[3], 1, mp ]
    for id in 1:dims[3]
        # cs = [ minColData, maxColData, minColUVW, maxColUVW ]
        cs = [ 1, IA[2], 1, mp ]
        for iw in 1:dims[2]
            # rs = [ minRowData, maxRowData, minRowUVW, maxRowUVW ]
            rs = [ 1, IA[1], 1, mp ]
            for ih in 1:dims[1]

                grid[count] = ( rs[1], rs[2], cs[1], cs[2], ds[1], ds[2],
                                rs[3], rs[4], cs[3], cs[4], ds[3], ds[4] )
                count += 1;
                @simd for idx in 1:4 ( rs[idx] += rdelta[idx] ) end
            end
            @simd for idx in 1:4 ( cs[idx] += cdelta[idx] ) end
        end
        @simd for idx in 1:4 ( ds[idx] += ddelta[idx] ) end
    end
    end #@inbounds

    return grid
end

# compute coordinates of search area after shifting them by previously computed displacements during multi-pass
function searchGrid( interGrid, SM, offh, offw, offd, h, w, d )
    r1 = max( interGrid[1] - SM[1] + offh, 1 );
    r2 = min( interGrid[2] + SM[1] + offh, h );
    c1 = max( interGrid[3] - SM[2] + offw, 1 );
    c2 = min( interGrid[4] + SM[2] + offw, w );
    d1 = max( interGrid[5] - SM[2] + offd, 1 );
    d2 = min( interGrid[6] + SM[2] + offd, d );
    so = max( 0, r1 - ( interGrid[1] - SM[1] ) - offh ),
         max( 0, c1 - ( interGrid[3] - SM[2] ) - offw ),
         max( 0, d1 - ( interGrid[5] - SM[3] ) - offd );

    return ( r1, r2, c1, c2, d1, d2 ), so
end



function setTo0!( A::Array{T,N} ) where {T<:Number,N}
    zer = T(0)
    @inbounds @simd for i in 1:length(A)
        A[ i ] = zer
    end
end

function setTo0!( A::Array{T,N}... ) where {T<:Number,N}
    setTo0!.( A )
end




# optimized to use SIMD, performing multiple comparisons at the same time
function maxval( a::AbstractArray{T,N} ) where {T<:Real,N}
    m1 = a[1]
    m2 = a[2]
    m3 = a[3]
    m4 = a[4]
    m5 = a[end]

    @inbounds begin
        @simd for i in 5:4:length(a)-4
            m1 = ( a[ i ] > m1 ) ? a[ i ] : m1
            m2 = ( a[i+1] > m2 ) ? a[i+1] : m2
            m3 = ( a[i+2] > m3 ) ? a[i+2] : m3
            m4 = ( a[i+3] > m4 ) ? a[i+3] : m4
        end
        for i in length(a)-4:length(a)-3
            m5 = ( a[i] > m5 ) ? a[i] : m5
        end
    end
    return maximum( [ m1, m2, m3, m4, m5 ] )
end




# optimized to perform parallel computations with simd
function firstPeak( A::Array{T,N} ) where {T<:Real,N}

    if length(A) < 12
        mv, mx = findmax( A )
        return Tuple(mx), mv
    end

    # maximum indices for parallel SIMD processing
    mi1 = 1
    mi2 = 2
    mi3 = 3
    mi4 = 4
    mi5 = length(A)

    @inbounds begin
    @simd for i in 5:4:length(A)-4
        mi1 = ( A[ i ] > A[mi1] ) ?  i  : mi1;
        mi2 = ( A[i+1] > A[mi2] ) ? i+1 : mi2;
        mi3 = ( A[i+2] > A[mi3] ) ? i+2 : mi3;
        mi4 = ( A[i+3] > A[mi4] ) ? i+3 : mi4;
    end

    for i in length(A)-4:length(A)-1
        mi5 = ( A[i] > A[mi5] ) ? i : mi5;
    end
    end #inbounds

    # Finding the max value among the "independently"/SIMD computed maximum indices
    val, idx = findmax( [  A[ mi1 ], A[ mi2 ], A[ mi3], A[ mi4 ], A[ mi5 ] ] );

    maxindex = (idx==1)*mi1 + (idx==2)*mi2 + (idx==3)*mi3 +
               (idx==4)*mi4 + (idx==5)*mi5;

    # Transforming from linear indexing to cartesian
    cidx = [ 0 for i in 1:N ];
    if N == 3
        h, w, d = size(A)
        z = ceil( Int64, maxindex/(h*w) )
        x = ceil( Int64, (maxindex - (z-1)*h*w)/h )
        y = maxindex - (x-1)*h - (z-1)*h*w;
        cidx = ( y, x, z );
    else
        h, w = size(A)
        x = ceil( Int64, maxindex/h )
        y = maxindex - (x-1)*h;
        cidx = ( y, x );
    end

    return cidx, A[ maxindex ]
end




function secondPeak( Arr::A{T,N}, rad::I, idx::NTuple{N,I}) where {T<:Real,N}

    # extract region around first peak, before setting it to 0
    mins   = max.(     1    , idx .- rad );
    maxs   = min.( size(Arr), idx .+ rad );
    ranges = [ mins[e]:maxs[e] for e in 1:N ];
    OGvals = copy( Arr[ ranges... ] );

    # setting region to 0
    Arr[ ranges... ] .= T(0);

    # computing first peak
    p2, val = firstPeak( Arr );

    # setting region back to its original value
    Arr[ ranges... ] .= OGvals

    return p2, val
end



# Used in subpixel approximation functions
function indexOnBorder( index::NTuple{N,I}, csize::NTuple{N,I} ) where {N}
    return any( index .== 1 ) || any( index .== csize )
end

neighIndices( idx::II  ) = [(idx[1]+y,idx[2]+x)   for y in -1:1, x in -1:1 ]
neighIndices( idx::III ) = [(idx[1]+y,idx[2]+x,idx[3]+z) for y in -1:1, x in -1:1, z in -1:1 ]

function getNeighbours( index::NTuple{N,I}, cmatrix::Array{T,N} ) where {T<:AbstractFloat,N}
    indices = neighIndices( index );
    neighs  = [ cmatrix[ x... ] for x in indices ]
    return neighs
end




function putWithinSearch!( search::Array{T,N}, A::Array{U,N},
                           so::NTuple{3,Int64}, coord::NTuple{6,Int64} 
                         ) where {T<:Real,U<:Real,N}
    start  = so .+ 1;
    r1, r2 = (coord[1], coord[2])
    c1, c2 = (coord[3], coord[4])
    d1, d2 = (N == 2) ? (1,1) : (coord[5], coord[6])

    @inbounds begin
    for z in 0:(d2-d1)
        for c in 0:(c2-c1)
            @simd for r in 0:(r2-r1)
                search[start[1]+r, start[2]+c, start[3]+z] = A[r1+r, c1+c, d1+z]
            end
        end
    end
    end # inbounds
end

function putWithinPadded!( padA::Array{Complex{T},N}, A::AbstractArray{U,N}, mean,
                           so::NTuple{3,Int64}, coord::NTuple{6,Int64}
                         )  where {T<:AbstractFloat,U<:Real,N}
    off    = so .+ 1;
    r1, r2 = (coord[1], coord[2])
    c1, c2 = (coord[3], coord[4])
    d1, d2 = (N == 2) ? (1,1) : (coord[5], coord[6])

    Tzero = T(0)
    Tmean = convert( T, mean )

    @inbounds begin
    for z in 0:(d2-d1)
        for c in 0:(c2-c1)
            @simd for r in 0:(r2-r1)
                padA[ off[1]+r,
                      off[2]+c,
                      off[3]+z ] = complex( convert( T, A[r1+r,c1+c,d1+z] ) - Tmean, Tzero )
            end
        end
    end
    end # inbounds
end

function putWithinPadded!( padA::Array{Complex{T},N}, A::AbstractArray{U,N}, mean 
                         )  where {T<:AbstractFloat,U<:Real,N}
    Tzero = T(0)
    Tmean = convert( T, mean )

    @inbounds begin
    for zet in 1:size( A, 3 )
        for col in 1:size( A, 2 )
            @simd for row in 1:size( A, 1 )
                padA[row,col,zet] = complex( convert( T, A[row,col,zet] ) - Tmean, Tzero )
            end
        end
    end
    end #@inbounds
end
