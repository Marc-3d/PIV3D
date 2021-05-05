# PIV3D
Free and open-source three-dimensional PIV implementation for research. 


### Installation instructions.

To install this package in Julia you should first install the LIBTIFF package, and then PIV3D: 
```Julia
using Pkg
Pkg.add(url="https://github.com/Marc-3d/LIBTIFF.jl.git")
Pkg.add(url="https://github.com/Marc-3d/PIV3D.git")
```

### How to run PIV3D

After installing the package, running a PIV analysis takes three steps.   

**1st** Load your data. (Look for packages that can open your data, ex [HDF5](https://github.com/JuliaIO/HDF5.jl)) <br>
**2nd** Create an object of custom PIV parameters. <br>
**3rd** Call the PIV function on two consecutive images/volumes.

```Julia
using PIV3D

# Thus far I have worked with volumes stored as .tiff files. I have modified LIBTIFF to read 3D scanline volumes. 
volume1 = PIV3D.PIVload( "filename1.tiff", typ=Float32, path="/path/to/data/" )
volume2 = PIV3D.PIVload( "filename2.tiff", typ=Float32, path="/path/to/data/" )

# interArea, overlap and searchMargin can be given different values for each dimension, ex interArea=(10,30,20), overlap=(5,5,8)...
params = PIV3D.setPIVParameters( interArea=32, overlap=10, searchMargin=5, corr="ZNCC", mpass=2 ); 
 
U, V, W, SN = PIV3D.PIV( volume1, volume2, params );
```

The most meaningful parameters to change are the interrogation area size (interArea), the overlap between adjacent interrogation areas (overlap), the search margin (searchMargin), the cross-correlation algorithm (corr) and multi-pass depth (mpass). The full list of options can be found in [the source code](src/parameters.jl). Feel free to contact me if you have any questions or suggestions about additional parameters that should be added to PIV3D.

### Visualizing the results

To visualize my 3D volumes and the 3D vector fields I have been using Paraview. For this reason, I have added functions to write these into VTK files. 
```Julia
vectorFieldToVTK( "filename", U, V, W, path="/path/to/your/results/" )
volumeToVTK( "filename", volume1, path="/path/to/your/liking/" )
```

### Evaluation instructions. 

The evaluation code is included in the Evaluation branch. 

1-. Clone the evaluation branch in your computer: 

```
git clone -b Evaluation https://github.com/Marc-3d/PIV3D.git
```

2-. Open a Julia terminal and install IJulia: 

```Julia
using Pkg
Pkg.add("IJulia")
using IJulia
```

NOTE: After running "using IJulia" for the first time, you will be asked whether IJulia should download and install miniconda. This is needed to run the jupyter notebook. In windows, I got this [error] which is solved by temporarily stopping Kasperky antivirus while miniconda is downloaded and installed. 

[error]: https://discourse.julialang.org/t/problem-with-curl-exe-windows-and-package-installation/29525/21

3-. Run the jupyter notebook from Julia: 

```Julia
IJulia.notebook()
```

4-. A tab will open in your default web browser. Make your way in the jupyter tree to the PIV3D folder that was downloaded in step 1. The "EvaluationNotebooks"  folder contains the Julia evaluation notebook, python performance evaluation notebook and .cpp performance script. 
