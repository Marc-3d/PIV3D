export evaluate, readResults, saveResults


function parametersToString( synParams::SP, pivParams::PP, tfmParams::TP )
    
    name = ""; 
    for field in fieldnames( SP )
        name = string( name, field, ":", getproperty( synParams, field ), ">" )
    end
    for field in fieldnames( PP )
        name = string( name, field, ":", getproperty( pivParams, field ), ">" )
    end
    for field in fieldnames( TP )
        name = string( name, field, ":", getproperty( tfmParams, field ), ">" )
    end
    return name[1:end-1]
end

function metaParametersToString( metaParams::MP )
    
    name = ""; 
    for field in fieldnames( metaParameters )
        name = string( name, field, ":", getproperty( metaParams, field ), ">" )
    end
    return replace( name[1:end-1], " " => "" )
end

function saveResults( synParams::SP, pivParams::PP, tfmParams::TP, metaParams::MP, results; savedir=""::String, tag=""::String )
    
    id       = parametersToString( synParams, pivParams, tfmParams )
    filename = metaParametersToString( metaParams )
    filepath = string( savedir, tag, filename, ".txt" )
    
    if ( isfile( filepath ) )
        io = open( filepath, "a")
         println( io, "" )
         println( io, id )
        writedlm( io, results )
        close(io)
    else
        io = open( filepath, "w")
         println( io, id )
        writedlm( io, results )
        close(io)
    end
end

function readResults( metaParams::MP; prefix="" )
    
    filename = metaParametersToString( metaParams ); 

    if ( isfile( string( prefix, filename, ".txt" ) ) )
        data = [ [], ]
        idx  = 1; 
        for line in eachline( string( prefix, filename, ".txt" ) )
            if ( occursin( ">", line ) )
                continue;
            elseif ( line == "" ) 
                idx += 1; 
                push!( data, []); 
                continue;
            else
                push!( data[idx], tuple( parse.(Float64, split( line, "\t") ) ... ) )
            end
        end
        return data
    end
    return false
end



"""
    evaluate( metaParams::MP, synParams::SP, pivParams::PP, tfmParams::TP )

    Calculates the biases and random errors of PIV vector fields calculated from 
    synthetic data as one parameters changes. A known translation is applied to the 
    synthetic particles, allowing to measure the mean difference (bias) between PIV 
    predictions and the ground truth, and the standard deviation of the PIV 
    predictions (random error). 

    metaParams: Determines the changing variable, its value interval and the number 
                of repeats each bias and random error calculation.
     
     synParams: Determines the amount and size of the synthetic Particles.
    
     pivParams: Determines the size of the search and interrogation area, the 
                overlap and mp_factor.
     
     tfmParams: Determines the transformation properties. Only translations are 
                supported by now. 
                Translations are sampled form a multivariate distribution with a 
                mean and variance and covariance. 

    All PIV evaluations shown in the paper can be found in "PIVevaluations.ipynb"
    jupyter notebook. 
"""
function evaluate( metaParams::MP, synParams::SP, pivParams::PP, tfmParams::TP; silent=false, checkParams=false )
    
    if synParams.d == 1
        tfmParams = to2DTransform( tfmParams ) # ensure tfmParams bears a 2D transform
    else 
        tfmParams = to3DTransform( tfmParams ) # ensure tfmParams bears a 3D transform 
    end
    return evaluate( metaParams, ( synParams, pivParams, tfmParams ), silent=silent, checkParams=checkParams ); 
end

function evaluate( metaParams::MP, parameters::Tuple{SP,PP,TP2D}; silent=false, checkParams=false )
    
    println("\n2D PIV evaluation of $(metaParams.variable)"); 
   
    index = 0;
    for idx in 1:3
        metaParams.variable in fieldnames( typeof(parameters[idx]) ) && ( index = idx; break; )
    end
    ( index == 0 ) && return results;
    
    og_val = getproperty( parameters[index], metaParams.variable ); # original value of changing variable
    valtype = ( typeof(og_val) <: Tuple ) ? eltype( og_val ) : typeof(og_val)
    
    # Array-like results -> results[i] = ( x_bias, y_bias, rand_ex, rand_ey ) i
    results = Array{NTuple{4,Float64},1}(undef,0); 
    
    # DataFrame results
    resultsDF = DataFrame(  
       var     = Array{ valtype, 1}(undef, length(collect(0:metaParams.metaloop))),
       x_bias  = collect( Float64, 0:metaParams.metaloop ), 
       y_bias  = collect( Float64, 0:metaParams.metaloop ), 
       x_error = collect( Float64, 0:metaParams.metaloop ), 
       y_error = collect( Float64, 0:metaParams.metaloop ) 
    );
    
    step   = ( metaParams.metaloop > 0 ) ? (metaParams.max - metaParams.min)/metaParams.metaloop : 0;
        
    for m in 0:metaParams.metaloop # OUTER LOOP: metaParams.variable is updated.
        
        new_val = one.( og_val ) .* metaParams.min .+ m*step;
        new_val = ( typeof( og_val ) <: Integer ) ? round( Int64, new_val ) : new_val; 
        
        updateParameters!( parameters[index], metaParams.variable, new_val ); 
                
        !silent && println("metaloop $m, $(typeof(parameters[index])).$(metaParams.variable) = $new_val")
         silent && print("$m,") 

        tmp_res = [ 0.0, 0.0, 0.0, 0.0 ]; # [ bias_x, bias_y, x_error, y_error ]
        meanUV  = [ 0.0, 0.0 ];           # [ mean(u), mean(v) ]
        N       = 0; 
        
        while N < metaParams.repeats
            
            # creating synthetic images
            trnsform = transformation( parameters[3] );
            parts    = generateParticles( parameters[1], pivparams=parameters[2] ); 
            parts_t  = applyTransformation(  parts, trnsform  ); 
            image1   = renderParticles( parts  , parameters[1] );
            image2   = renderParticles( parts_t, parameters[1] );

            # out[1] = U, out[2] = V, out[3] = sigNoise
            out = PIV( image1, image2, parameters[2], checkParams=checkParams ); 
            
            sizeUV = size(out[1]);
            
            N += (sizeUV[2] - 2)*(sizeUV[1] - 2); # TODO: ensure N > 1
            
            meanUV[1] = Statistics.mean(out[1]);
            meanUV[2] = Statistics.mean(out[2]); 
            
            for i in 2:(sizeUV[2] - 1)
                for j in 2:(sizeUV[1] - 1)
                    for e in 1:2
                        tmp_res[ e ] += abs(out[e][j,i] - trnsform[e]); # biases
                        tmp_res[e+2] += ( meanUV[e] - out[e][j, i] )^2; # rand_err
            end end end
        end
        @inbounds for e in 1:4 ( tmp_res[e] /= N ) end
        
        resultsDF[ m+1, :var     ] = ( typeof(new_val) <: Tuple ) ? new_val[1] : new_val;
        resultsDF[ m+1, :x_bias  ] = tmp_res[1]; 
        resultsDF[ m+1, :y_bias  ] = tmp_res[2]; 
        resultsDF[ m+1, :x_error ] = sqrt(tmp_res[3]);
        resultsDF[ m+1, :y_error ] = sqrt(tmp_res[4]); 
        
        push!( results, (    tmp_res[1] ,      tmp_res[2], 
                        sqrt(tmp_res[3]), sqrt(tmp_res[4]) ) )
    end
    
    # setting variable back to OG_val
    setproperty!( parameters[index], metaParams.variable, og_val ); 
    
    println( "\nvariable values: ", [ round( x[1], digits=3 ) for x in resultsDF.var ] ); 
    
    return results, resultsDF
end

function evaluate( metaParams::MP, parameters::Tuple{SP,PP,TP3D}; silent=false, checkParams=false )
    
    println("\n3D PIV evaluation of $(metaParams.variable)"); 
    
    index = 0;
    for idx in 1:3
        if metaParams.variable in fieldnames( typeof(parameters[idx]) )
            index = idx; break;
        end
    end
    ( index == 0 ) && return results;
    
    og_val = getproperty( parameters[index], metaParams.variable );
    valtype = ( typeof(og_val) <: Tuple ) ? eltype( og_val ) : typeof(og_val)
    
    # Arary-like results: 
    # results[i] = ( x_bias, y_bias, z_bias, rand_ex, rand_ey, rand_ez )
    results = Array{NTuple{6,Float64},1}(undef,0); 
    
    # DataFrame results
    outerloop = 0:metaParams.metaloop; 
    outerlen  = length(collect(outerloop));
    
    resultsDF = DataFrame(  
       var     = Array{ valtype, 1 }( undef, outerlen ),
       x_bias  = collect( Float64, outerloop ), 
       y_bias  = collect( Float64, outerloop ),
       z_bias  = collect( Float64, outerloop ),
       x_error = collect( Float64, outerloop ), 
       y_error = collect( Float64, outerloop ),
       z_error = collect( Float64, outerloop )
    );
                    
    # quantity by which the changing variable is updated each iteration
    step = ( metaParams.max - metaParams.min ) / metaParams.metaloop; 
    ( metaParams.metaloop <= 0 ) && ( step = 0; )
    
    for m in outerloop
        
        new_val = one.( og_val ) .* metaParams.min .+ m * step;
        new_val = (typeof(og_val)<:Integer) ? round(Int64,new_val) : new_val; 
        
        updateParameters!( parameters[index], metaParams.variable, new_val ); 
        
        if silent
            print("$m,") 
        else
            println("metaloop $m, 
            $(typeof(parameters[index])).$(metaParams.variable) = $new_val")
        end

        tmp_res  = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ];
        meansUVW = [ 0.0, 0.0, 0.0 ];      
        N        = 0; 

        while N < metaParams.repeats
                        
            trnsform =    transformation( parameters[3] ); 
            parts    = generateParticles( parameters[1], pivparams=parameters[2] ); 
            parts_t  = applyTransformation( parts,   trnsform    ); 
            volume1  = renderParticles( parts    , parameters[1] );
            volume2  = renderParticles( parts_t  , parameters[1] );
                     
            # out[1] = U, out[2] = V, out[3] = W, out[4] = SN
            out = PIV( volume1, volume2, parameters[2], checkParams=checkParams ); 
            
            sizeUVW = size(out[1]); 
            
            N += (sizeUVW[1] - 2)*(sizeUVW[2] - 2)*(sizeUVW[3] - 2)
            
            meansUVW[1] = Statistics.mean(out[1]); 
            meansUVW[2] = Statistics.mean(out[2]); 
            meansUVW[3] = Statistics.mean(out[3]);
            
            for k in 2:(sizeUVW[3] - 1)
                for i in 2:(sizeUVW[2] - 1)
                    for j in 2:(sizeUVW[1] - 1)
                        for e in 1:3
                            tmp_res[ e ] += abs(out[e][j,i,k] - trnsform[e]); 
                            tmp_res[e+3] += (meansUVW[e] - out[e][j,i,k])^2;
            end end end end
        end
        @inbounds for e in 1:6 ( tmp_res[e] /= N ) end
            
        resultsDF[ m+1, :var     ] = ( typeof(new_val) <: Tuple ) ? new_val[1] : new_val; 
        resultsDF[ m+1, :x_bias  ] = tmp_res[1]; 
        resultsDF[ m+1, :y_bias  ] = tmp_res[2]; 
        resultsDF[ m+1, :z_bias  ] = tmp_res[3];
        resultsDF[ m+1, :x_error ] = sqrt(tmp_res[4]);
        resultsDF[ m+1, :y_error ] = sqrt(tmp_res[5]); 
        resultsDF[ m+1, :z_error ] = sqrt(tmp_res[6]); 
        
        push!( results, (    tmp_res[1] ,      tmp_res[2] ,      tmp_res[3], 
                        sqrt(tmp_res[4]), sqrt(tmp_res[5]), sqrt(tmp_res[6]) ) )
    end
    
    # setting variable back to OG_val
    updateParameters!( parameters[index], metaParams.variable, og_val ); 
    
    println( "\nvariable values: ", [ round( x[1], digits=3 ) for x in resultsDF.var ] ); 
    
    return results, resultsDF
end

using Measures

"""
    This function contains the styling used to generate the figures in the paper. 
    Plots are produced using Gadfly.
"""
function prepareLayer(  field::Symbol, _min, _max, 
                        results::Array{NTuple{4,Float64},1};  
                        fun=(x->x), color="blue" ) 
    
    index = findmax( [ field == :bias_x , field == :bias_y, 
                       field == :rand_ex, field == :rand_ey, true ] )[2];
    
    ( index == 5 ) && error("pass one of :bias_x, :bias_y, :rand_ex or :rand_ey");
    
    y = [ fun( x[index] ) for x in results  ];
    x = _min .+ collect( 0:length( results )-1 )./length(results) .* _max; 
            
    return Gadfly.layer( x=x, y=y, Geom.line, Theme(default_color=color,) ); 
end

function plotResults(  field::Symbol, _min, _max, 
                       results::Array{NTuple{4,Float64},1}...; 
                       fun=(x->x), xlabel="", ylabel="" ) 
    
    
    colors = [ "red", "cyan3", "blue", "orange", "gold" ]
    
    layers = Array{Gadfly.Layer,1}(undef, length(results) );
    
    for i in 1:length(results)
        
       layers[i] = prepareLayer( field, _min, _max, results[i], 
                                 fun=fun, color=colors[i] )[1]; 
    end
    
    p1 = Gadfly.plot( Coord.Cartesian( xmin=_min, xmax=_max), 
                      Guide.xlabel(xlabel),
                      Guide.ylabel(ylabel, :vertical),
                      layers... );  
end    
    
