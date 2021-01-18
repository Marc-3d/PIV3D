export plotVolume

function plotVolume( vol; gui=false)
    
    # If the axes don't have the same size, PyPlot freaks out. Therefore, I get the largest dimensions, and make 
    # a cube out of it. Then, volume will be copied inside the cube, normalizing values between 0 and 1.
    
    larg = max( size(vol,1), size(vol,2), size(vol,3) ); 
    cube = zeros( Float32, larg, larg, larg ); 
    
    MIN, MAX = minimum(vol), maximum(vol); 
    for z in 1:size(vol,3)
        for x in 1:size(vol,2)
            for y in 1:size(vol,1)
                cube[y,x,z] = (vol[y,x,z] - MIN )/(MAX-MIN)
            end
        end
    end
    
    color = [ (0,0,0,x) for x in cube  ];
    
    # ADDING BORDERS AROUND the original volume
    for z in 1:size(vol,3)
        color[ size(vol,1), size(vol,2), z ] = ( 1, 0, 0, 0.5 ); 
        color[ size(vol,1),           1, z ] = ( 1, 0, 0, 0.5 ); 
        color[           1, size(vol,2), z ] = ( 1, 0, 0, 0.5 ); 
        color[           1,           1, z ] = ( 1, 0, 0, 0.5 ); 
    end
    for x in 1:size(vol,1)
        color[ x,           1, size(vol,3) ] = ( 1, 0, 0, 0.5 ); 
        color[ x, size(vol,2), size(vol,3) ] = ( 1, 0, 0, 0.5 ); 
        color[ x,           1,           1 ] = ( 1, 0, 0, 0.5 ); 
        color[ x, size(vol,2),           1 ] = ( 1, 0, 0, 0.5 ); 
    end
    for y in 1:size(vol,2)
        color[           1, y, size(vol,3) ] = ( 1, 0, 0, 0.5 ); 
        color[ size(vol,1), y, size(vol,3) ] = ( 1, 0, 0, 0.5 ); 
        color[           1, y,           1 ] = ( 1, 0, 0, 0.5 ); 
        color[ size(vol,1), y,           1 ] = ( 1, 0, 0, 0.5 ); 
    end
    
    xs  = [ x for x in 1:size(cube,1), y in 1:size(cube,2), z in 1:size(cube,3) ]
    ys  = [ y for x in 1:size(cube,1), y in 1:size(cube,2), z in 1:size(cube,3) ]
    zs  = [ z for x in 1:size(cube,1), y in 1:size(cube,2), z in 1:size(cube,3) ]
    
    f = figure(); 
    PyPlot.pygui( gui ); 
    PyPlot.scatter3D( xs, ys, zs, c=color[:] );   
end

function volumeToVTK( vol; filename="vol.vtk" )
    
    w, h, d = size(vol); 
    
    txt="""# vtk DataFile Version 4.2
    Vtk file written from PIV3D.jl
    ASCII
    DATASET STRUCTURED_POINTS
    DIMENSIONS $w $h $d
    ORIGIN 0 0 0
    SPACING 1 1 1
    
    POINT_DATA $(length(vol))
    SCALARS intensity double 1
    LOOKUP_TABLE default""" 
        
    for e in vol 
        txt *= "\n$e"
    end
    
    io = open( filename, "w")
    println( io, txt )
    close(io)
end