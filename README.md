# Triple Gravitational Microlens Library
This repository holds a library consisting of various modelling tools for [gravitational microlensing](https://en.wikipedia.org/wiki/Gravitational_microlensing). Most of the ideas included are products of my PhD research and this project aims to both put the algorithms involved to one place and transcribe them into tidy readable/maintainable object-oriented c++ code. 

## Installation Notes
Running 'make' in the root directory should compile the source code a produce binary 'ccc\_test' and shared library 'ccc.so'. The compilation was tested on Fedora 29 and Ubuntu 20 distributions.

## Main Components of the library

### Laguerre's Method
My code uses polynomial solving extensively. For that I use [Laguerre's Method](https://en.wikipedia.org/wiki/Laguerre%27s_method). My version of the algorithm is greatly inspired by book 'Numerical Recipes 3rd Edition: The Art of Scientific Computing', however, it is implemented with emphasis on re-using solution from previous step in sequence of changes for a parametric polynomial, e.g., Critical Curves and Light Curves.

### CCC 
Class for producing Critical Curves and Caustics, which I use to characterise the microlenses. The class also includes a wrapper to be used within a python.

### ImgPoint 
Class for simulating images produced by triple-lens microlensing event, also comes with python wrapper.

### LightCurves
Atributes amplification to images and outputs them as a time sequence for corresponding source position. In basic form it counts amplification of point sources which is sum inverse values of the lens-equation Jacobian in the image position for each image position. 
The derived class, LightCurveIRS, where IRS stands for 'inverse ray shooting' simulated images produced by extended-source triple-lens microlensing. For that purpose a flood filling algorithm is employed so that each image is filled to its full size determined by shooting rays backward from the image to source.
LightCurveIRS has a python wrapper.  

### Examples
Contains python scripts that call shared library files to execute compiled c++ code. Script 'run\_cpplib.py' to call critical-curve functionality, 'run\_imglib.py' to call ImgPoint, 'run\_lclib.py' to call LightCurves.

### main.cc
As this is mostly meant as a shared library to be later used within a Python application, main.cc serves for regression testing. Plan is to move to unit testing as the library grows larger. 

## Unit tests
The project utilizes [Google Tests](https://github.com/google/googletest/) for unit testing.
To produce the test binaries, download the content of the [test repository](https://github.com/google/googletest/) into '/test/googletest' directory and run CMake, e.g., by running 'cmake .. -DCMAKE\_BUILD\_TYPE=Release -G \"Unix Makefiles\"' in the '/bin' directory. This will produce a Makefile, to produce a binaries run 'make' in the '/bin' directory. Resulting test binary will be made in 'bin/test' directory.

If you don't have CMake installed run, e.g., 'sudo yum cmake' (Fedora), 'sudo apt-get install cmake' (Ubuntu).

## CUDA functionality

In '/src/Cuda' directory is functionality utilizing CUDA for inverse-ray shooting on GPU. To build the code to GPU, use 'make cuda'. Alternatively, if you want to make sure that the GPU-capable code won't get included, use 'make cpu'.


## Tools worth mentioning

### GPROF

A simple profiling tool that gives a good overview of time spent on each function. 
Make an executable by running 'make profile', then run the resulting executable in the 'bin/' directory. Finally, run 'gprof ccc\_profile gmon.out'. 
