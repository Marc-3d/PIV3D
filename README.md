# PIV3D
Free and open-source three-dimensional PIV implementation for research. 

### Installation instructions.

To install this package in Julia you should first install the LIBTIFF package, and then PIV3D: 
```Julia
using Pkg
Pkg.add(url="https://github.com/Marc-3d/LIBTIFF.jl.git")
Pkg.add(url="https://github.com/Marc-3d/PIV3D.git")
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
