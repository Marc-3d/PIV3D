include("znccUtils.jl")

#= 
	2D IMPLEMENTATION 
=#


#  Lewis, J.P.: Fast normalized cross-correlation. Ind. Light Magic10(2001)
function crossCorrelation!( ::ZNCC, ccr::A{T,2}, num::A{T,2},
                                    padF::A{C{T},2}, padG::A{C{T},2}, plan, iplan,
                                    sumF::A{U,2}, sumF2::A{U,2}, meanF, meanF2, 
                                    stdG, sizeF, sizeG, shifts
                          ) where {T<:Real,U<:Real}
    
    # Computing numerator with FFT.
    crossCorrelation!( FFT(), padF, padG, ccr, plan, iplan );
    Base.circshift!( num, ccr, shifts ) 
    
    # Computing denominator through summed-area tables.
    N = prod(sizeG); 
    N2F, stdI = N*meanF2, stdG*sqrt(N)
    @inbounds for col in 1:size(ccr,2)
        for row in 1:size(ccr,1)
            F, F2 = znccSumValues( sumF, sumF2, row, col, sizeF, sizeG )
            ccr[ row, col ] = num[ row, col ]/( sqrt( F2 + N2F - 2*F*meanF )*stdG );
        end
    end
end

function PIV_2D( ::ZNCC, img1::A{T,2}, img2::A{T,2},
                         IA::III, SM::III, overlap::III, mpass::I, width::I,
                         peak::S, sigNoise::S, filtFun::Function, threshold::F, 
                         mpFun::Function, reps::Array{I,1}; 
                         corrType=Float32 
			   ) where {T<:Real}  
    
	# Calculating vector field size.
    h, w   = size( img1 ); 
    step   = IA .- overlap; 
    VFsize = 1 .+ floor.( Int64, ((h,w,0).-IA)./step ); 
    U  = zeros( corrType, VFsize[1:2] ); 
    V  = zeros( corrType, VFsize[1:2] ); 
    SN = zeros( corrType, VFsize[1:2] );

	ignoreSN = sigNoise == ""; 
    
    @inbounds for mp in mpass:-1:1
                
        # Scaling IA, SM and step (overlap) to multi-pass iteration.
        IA_mp   = mpFun.(  IA , mp ); 
        SM_mp   = mpFun.(  SM , mp ); 
        step_mp = mpFun.( step, mp ); 
        last_mp = ( mp == 1 ); 
   
        ssize   = IA_mp .+ 2 .* SM_mp
        search  = zeros( T, ssize[1:2] ); 
        sumS    = zeros( T, ssize[1:2] .+ 1 ); 
        sumS2   = zeros( T, ssize[1:2] .+ 1 ); 

        csize   = 2 .* ( IA_mp .+ SM_mp  );
        cmatrix = zeros( corrType, csize[1:2] );
        num     = zeros( corrType, csize[1:2] ); 
        pads    = zeros( Complex{corrType}, csize[1:2] ); 
        padi    = zeros( Complex{corrType}, csize[1:2] ); 
        plan    =  FFTW.plan_fft!( padi ); 
        iplan   = FFTW.plan_ifft!( padi );
		shifts  = div.( csize, 2 )[1:2] .+ SM_mp[1:2] .- 1;
            
        # Constructing grid coordinates for each interrogation volume, 
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

            # 1-. Shifting the searchArea ranges by previously calculated displacements 
            offH = round(Int64, U[ g[7], g[9] ]);
            offW = round(Int64, V[ g[7], g[9] ]);
            
            # 1.2-. Search volume coordinates after shifting.   
            scoords, so = searchGrid( g, SM, offH, offW, 0, h, w, 0 )

            setTo0!( search )
            putWithinSearch!( search, img2, so, scoords )

            setTo0!( sumS  ); 
            setTo0!( sumS2 ); 
            znccSumArrays!( search, sumS, sumS2 ); 

            meanI = Statistics.mean( interr );
            meanS = Statistics.mean( search );

            setTo0!( pads )
            setTo0!( padi )
            putWithinPadded!( padi, interr, meanI ); 
            putWithinPadded!( pads, search, meanS );

            meanS2 = meanS*meanS; 
            meanI2 = meanI*meanI; 
            stdI   = Statistics.std( interr, mean=meanI );

            # 2-. Cross-Correlation
            crossCorrelation!( ZNCC(), cmatrix, num, pads, padi, plan, iplan, 
                               sumS, sumS2, meanS, meanS2, stdI, ssize, IA_mp, shifts );

            # 3-. Calculation of displacement
            ( r, c ) = approxTranslation( cmatrix, peak, Val(last_mp) )
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

function crossCorrelation!( ::ZNCC, ccr::A{T,3}, num::A{T,3},
                                    padF::A{C{T},3}, padG::A{C{T},3}, plan, iplan,
                                    sumF::A{U,3}, sumF2::A{U,3}, meanF, meanF2,
                                    stdG, sizeF, sizeG, shifts 
						  ) where {T<:Real, U<:Real}
    
    # Computing numerator with FFT.
    crossCorrelation!( FFT(), padF, padG, ccr, plan, iplan );  
    Base.circshift!( num, ccr, shifts )

    # Computing denominator through summed-area tables.
    N = prod(sizeG);
    N2F, stdG  = N*meanF2, stdG*sqrt(N);
    @inbounds for z in 1:size(ccr,3)
        for c in 1:size(ccr,2)   
            for r in 1:size(ccr,1)
                F, F2 = znccSumValues( sumF, sumF2, r, c, z, sizeF, sizeG )
                ccr[r,c,z] = num[r,c,z]/( sqrt( F2 + N2F - 2*F*meanF )*stdG );
            end
        end
    end
end

function PIV_3D( ::ZNCC, vol1::A{T,3}, vol2::A{T,3},
                         IA::III, SM::III, overlap::III, mpass::I, width::I,
                         peak::S, sigNoise::S, filtFun::Function, threshold::F, 
                         mpFun::Function, reps::Array{I,1}; 
                         corrType=Float32 
			   ) where {T<:Real} 
    
	# Calculating vector field size.
    h, w, d = size( vol1 ); 
    step   = IA .- overlap; 
    VFsize = 1 .+ floor.(Int64,(size(vol1).-IA)./step); 
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
   
        ssize   = IA_mp .+ 2 .* SM_mp
        search  = zeros( T, ssize ); 
        sumS    = zeros( T, ssize .+ 1 ); 
        sumS2   = zeros( T, ssize .+ 1 ); 

        # Initialize cross-correlation variables once per multi-pass iteration
        csize   = 2 .* ( IA_mp .+ SM_mp  );
        cmatrix = zeros( corrType, csize );
        num     = zeros( corrType, csize ); 
        pads    = zeros( Complex{corrType}, csize ); 
        padi    = zeros( Complex{corrType}, csize ); 
        plan    =  FFTW.plan_fft!( pads ); 
        iplan   = FFTW.plan_ifft!( pads );
	    shifts  = div.( csize, 2 ) .+ SM_mp .- 1;
            
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

            if threshold > 0 && filtFun( interr ) < threshold 
                continue;
            end

            # 1-. Shifting the searchArea ranges by previously calculated displacements
            offH = round(Int64, U[ g[7], g[9], g[11] ]);
            offW = round(Int64, V[ g[7], g[9], g[11] ]);
            offD = round(Int64, W[ g[7], g[9], g[11] ]);
                
            # 1.2-. Search volume coordinates after shifting.
            scoords, so = searchGrid( g, SM, offH, offW, offD, h, w, d )

            setTo0!( search )
            putWithinSearch!( search, vol2, so, scoords )

            setTo0!( sumS  ); 
            setTo0!( sumS2 ); 
            znccSumArrays!( search, sumS, sumS2 ); 

            meanI = Statistics.mean( interr );
            meanS = Statistics.mean( search );

            setTo0!( pads )
            setTo0!( padi )
            putWithinPadded!( padi, interr, meanI ); 
            putWithinPadded!( pads, search, meanS );

            meanS2 = meanS*meanS; 
            meanI2 = meanI*meanI; 
            stdI   = Statistics.std( interr, mean=meanI );

            # 2-. Cross-Correlation
            crossCorrelation!( ZNCC(), cmatrix, num, pads, padi, plan, iplan, 
                               sumS, sumS2, meanS, meanS2, stdI, ssize, IA_mp, shifts );

            # 3-. Calculation of displacement
            ( r, c, z ) = approxTranslation( cmatrix, peak, Val(last_mp) )
			if ( !ignoreSN && last_mp )
            	SN[ g[7], g[9], g[11] ] = SNRatio( cmatrix, sigNoise, width )
			end

            # 4-. Updating U, V, W matrices 
            t1 =  g[7]:g[8]; 
            t2 =  g[9]:g[10]; 
            t3 = g[11]:g[12]; 
            U[ t1, t2, t3 ] .= offH - r;
            V[ t1, t2, t3 ] .= offW - c;
            W[ t1, t2, t3 ] .= offD - z; 
        end
    end # Multi-pass loop
    
    return U, V, W, SN
end



""" ZNCC out of place CROSS CORRELATION """
function crossCorrelation( ::ZNCC, f::A{T,3}, g::A{T,3}; typ=Float32 ) where {T<:Real, U<:Real}
    
    sizef  = size( f ); 
    sizeg  = size( g ); 
    
    sumf   = zeros( typ, sizef .+ 1 ); 
    sumf2  = zeros( typ, sizef .+ 1 ); 
    znccSumArrays!( f, sumf, sumf2 );     
    
    csize = sizef .+ sizeg; 
    corr  = zeros( typ, csize ); 
    num   = zeros( typ, csize ); 
    padf  = zeros( Complex{typ}, csize ); 
    padg  = zeros( Complex{typ}, csize );
    plan  =  plan_fft!( padf ); 
    iplan = plan_ifft!( padf ); 
	shifts = div( csize,2 ) .- 1; 
    
    meanf = Statistics.mean( f ); 
    meang = Statistics.mean( g ); 
    
    putWithinPadded!( padf, f, meanf ); 
    putWithinPadded!( padg, g, meang );
    
    meanf2 = meanf^2; 
    meang2 = meang^2; 
    stdg   = Statistics.std( g, mean=meang );
    
    crossCorrelation!( ZNCC(), corr, num, padf, padg, plan, iplan, 
                       sumf, sumf2, meanf, meanf2, stdg, sizef, sizeg, shifts );
    
    GC.gc(); 
    return corr; 
end

function crossCorrelation( ::ZNCC, f::A{T,2}, g::A{T,2}; typ=Float32 ) where {T<:Real, U<:Real}
    
    sizef  = size( f ); 
    sizeg  = size( g ); 
    
    sumf   = zeros( typ, sizef .+ 1 ); 
    sumf2  = zeros( typ, sizef .+ 1 ); 
    znccSumArrays!( f, sumf, sumf2, (0,0) );     
    
    csize = sizef .+ sizeg; 
    corr  = zeros( typ, csize ); 
    num   = zeros( typ, csize ); 
    padf  = zeros( Complex{typ}, csize ); 
    padg  = zeros( Complex{typ}, csize );
    plan  =  plan_fft!( padf ); 
    iplan = plan_ifft!( padf ); 
	shifts = div( csize[1:2],2 ) .- 1; 
    
    meanf = Statistics.mean( f ); 
    meang = Statistics.mean( g ); 
    
    putWithinPadded!( padf, f, meanf ); 
    putWithinPadded!( padg, g, meang );
    
    meanf2 = meanf^2; 
    meang2 = meang^2; 
    stdg   = Statistics.std( g, mean=meang );
    
    crossCorrelation!( ZNCC(), corr, num, padf, padg, plan, iplan, 
                       sumf, sumf2, meanf, meanf2, stdg, sizef, sizeg, shifts);
    
    GC.gc(); 
    return corr; 
end


