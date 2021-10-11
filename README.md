# PFMICE

PFMICE
![alt text](https://github.com/zhangwq1990/PFMICE/blob/11c0ed452340c213a93eb67bba940eaed6a01564/logo.png)
## Description

This is a 3D test case of the code PFMICE

## Getting Started

### Dependencies

* PFMICE is written with Fortran. Before you try to install PFMICE on your local computer or
workstation, Fortran compiler (GCC or Intel Fortran) and MPI libraries (e.g., openmpi or Intel
MPI) must be installed. OpenMPI is freely available online:

https://www.open-mpi.org/software/ompi/v4.1/.

and it can be installed in the following steps:
```
gunzip -c openmpi-4.1.1.tar.gz | tar xf -
```
```
cd openmpi-4.1.1
```

```
./configure --prefix=/usr/local
```
```
<...lots of output...>
```
```
make all install
```

Intel oneAPI Toolkits can be found in:

https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html

### Installing

* A Makefile can be found in the folder. If your current compiling environment is GCC+OpenMPI,
in the Makefile, you should use:

```
FC = mpif90
```
```
FFLAGS := -O3 -ffixed-line-length-none -mcmodel=large -fdefault-real-8 -cpp -Wall -fcheck=all
```
For Intel compiler+Intel MPI, the setup in the Makefile will be:
```
FC = mpiifort
```
```
FFLAGS := -r8 -fpconstant -O3 -132 -cpp
```
To compile the code, you need to compile the libraries at first with:
```
make libraries
```
and then compile the source code with:
```
make
```
you can remove the executables with:
```
make clean
```



* Any modifications needed to be made to files/folders

### Executing program

* To run the code:
```
mpirun -np N PFMICE
```
where N is the number of CPU cores.

## Help

Any advise for common problems or issues.
You can drop an e-mail to zhangwq1990@yahoo.com

## Authors

Contributors names and contact info

Developer: Wenqiang Zhang

Email: zhangwq1990@yahoo.com

Co-Developer: Armin Shahmardi, Ziyang Huang, Xuerui Mao

## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the GUN general Licenses - see the LICENSE.md file for details


## References
<a id="1">[1]</a> 
Ziyang Huang, Guang Lin, and Arezoo M Ardekani. (2021). 
A consistent and conservative phase-field model for thermo-gas-liquid-solid flows including liquid-solid phase change.
arXiv preprint,arXiv:2102.06863, 2021.

<a id="1">[2]</a>
Michael S Dodd and Antonino Ferrante. 
