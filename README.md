## Coursework on numerical methods, 5 semester FAMI NSTU
Implementation of the finite element method (FEM) for non-stationary parabolic boundary value problem in two-dimensional space, cylindrical coordinates. Finite elements are triangles with linear basis functions. The time scheme is four-layer implicit.

## Introduction
Chosen differential equation system describes the multi-phase flow of an incompressible unstirred isothermal fluid in an open domain. The goal is to find an approximation of the pressure function. The function u can be approximated by a function $u_h$ using linear combinations of basis functions $u \approx u_h $  
$ u_h = \[ \sum_{i} u_i \psi_i \]$  
Here, $\psi_i$ denotes the basis functions and $u_i$ denotes the coefficients of the functions that approximate u with $u_h$.
As a result of this program, we will get the coefficients $u_i$.
# Installation
```
 git clone https://github.com/verybluebird/MFE_coursework
```
Excecute this project in MS Visual Studio.

## Modules description:
**Draw** - Phyton module for drawing a grid. It is needed to check the correctness of the grid generation.   
![grid](https://github.com/verybluebird/MFE_coursework/blob/main/grid.png?raw=true)
**grid_triangles_linear** - module for generating a grid for finite elements of the form: triangles with linear basic functions.  
**MFE** - module that takes grid from previous module and approximates the unknown function over the domain with given boundary conditions using Finite Element Method.       
# Layout
All modules  will have a similar directory and file layout

    main.cpp
    setup.h
    include/
      functions.h
      ...
    src/
      functions.cpp
      ...
    data/
      ...
    ...

The `include/` directory contains header files.

The `src/` directory contains implementations of the functions.

The `data/` directory contains input data for program.

## Testing
* Describe domain in data files of **MFE** module (see the description below). 
* To test  program on a specific desired function, specify a function as
   ```double u(double r, double z) { }``` in **grid_triangles_linear** module. This is needed to generate the first boundary conditions.
* Build **grid_triangles_linear** module and then run it in the `data/` directory of **MFE** module.
* To check the correctness of the grid generation run **Draw** module in the same directory.
* Set the coefficients of differential equation in **MFE** module.
 
* If desired function is given, specify this function as  
  ```double mfe::u_t(double r, double z, double t){ } ```  
  then calculate the right part of the differential equation and specify as
  ```double mfe::rightPart(int field, double r, double z, double t)	{}``` in **MFE** module.
* Build and run **MFE** module. As a result you will have `q.txt` file where function values at grid nodes in different time layers are stored. When programm is tested on spesific desired function The calculation error is also output to the file.



## Data input

**Draw** module  
Input files:
   * `elem.txt` - node coordinates
   
    r_coordinate z_coordinate
  The first node is the first line, the second node is the second, and so on.
   * `node.txt` - node numbers in the order of the elements
    node_1 node_2 node_3
  The first element is the first line, the second element is the second, and so on.  

**grid_triangles_linear** module  
Input files: 
* `domain_rz.txt` - description of the computational domain
 ```
  Rw - initial coordinate in r (left)  
  Rb - final coordinate in r (right)  
  NL - number of layers in z  
  H_i - thickness of layer i in z  
  K_i - structural permeability of layer i  
  Phi_i - layer i porosity  
  S2_i - oil saturation of layer i  
```

* `mesh_rz.txt` - grid settings

Let the z-coordinate grid be uniform in each layer (kz=1)  
``` 
  nr - number of partitions over r  
  nz - number of partitions over z  
  kr - discharging coefficient for r  
  M - grid nesting in space 
```

* `mesh_t.txt` - time grid settings  
```
  t1 - start time  
  t2 - end time  
  n_t - number of time layers  
  k_t - coefficient of discharging over time  
  M - grid nesting by time  
```
* `zone_perf_rz.txt` - description of the position and power of perforation zones
```
  Nzp - number of perforation zones  
  Pu - upper coordinate of the perforation zone  
  Pd - lower coordinate of the perforation zone  
  Tetta - power of the perforation zone  
```
* `phaseprop.txt` - phase properties
```
  Nph - number of phases  
  Mu - viscosity of phase i  
```
* `plast.txt` - reservoir pressure (on the border on the right)
```
  Pressure
```




