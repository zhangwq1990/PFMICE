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


## Demos

Solidification of a liquid (the density ratio of the liquid to the solid is 2:1):

![Alt Text](https://github.com/zhangwq1990/PFMICE/blob/581fee55f495da315892b3d34b0d5dc611d14cb3/demo/111.gif)

Water on a hydrophilic surface:
![Alt Text](https://github.com/zhangwq1990/PFMICE/blob/581fee55f495da315892b3d34b0d5dc611d14cb3/demo/333.gif)
Water on a hydrophobic surface:
![Alt Text](https://github.com/zhangwq1990/PFMICE/blob/581fee55f495da315892b3d34b0d5dc611d14cb3/demo/444.gif)

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
Shahmardi, Armin and Rosti, Marco Edoardo and Tammisola, Outi and Brandt, Luca,
A fully Eulerian hybrid immersed boundary-phase field model for contact line dynamics on complex geometries,
Journal of Computational Physics, 110468, 2021


<a id="2">[2]</a> 
Ziyang Huang, Guang Lin, and Arezoo M Ardekani. (2021). 
A consistent and conservative phase-field model for thermo-gas-liquid-solid flows including liquid-solid phase change.
arXiv preprint,arXiv:2102.06863, 2021.

<a id="3">[3]</a>
Michael S Dodd and Antonino Ferrante. 
A fast pressure-correction method for incompressible two-fluid flows. 
Journal of Computational Physics, 273:416–434, 2014.

<a id="4">[4]</a>
Huang, Z., Lin, G. and Ardekani, A.M., 2020. 
Consistent, essentially conservative and balanced-force phase-field method to model incompressible two-phase flows. 
Journal of Computational Physics, 406, p.109192.


<a id="5">[5]</a>
Jim Douglas. 
Alternating direction methods for three space variables. Numerische Mathematik,
4(1):41–63, 1962.

<a id="6">[6]</a>
N Li and S Laizet. 
2decomp&fft–a highly scalable 2d decomposition library and fft interface,
cray user group 2010 conference, edinburgh. 
URL http://www. 2decomp. org/pdf/17B-CUG2010-paper-Ning Li. pdf, 2010.

<a id="7">[7]</a>
Xianmin Xu, Yana Di, and Haijun Yu. 
Sharp-interface limits of a phase-field model with a generalized navier slip boundary condition for moving contact lines. 
Journal of Fluid Mechanics, 849:805–833, 2018.

<a id="8">[8]</a>
Ziyang Huang, Guang Lin, and Arezoo M Ardekani. 
Consistent and conservative scheme for incompressible two-phase flows using the conservative allen-cahn model. 
Journal of Computational Physics, 420:109718, 2020.
