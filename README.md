# Triple Gravitational Microlens Library
This repository holds a library consisting of various modelling tools for gravitational microlensing[gravitational microlensing](https://en.wikipedia.org/wiki/Gravitational_microlensing). Most of the ideas included are products of my PhD research and this project aims to both put the algorithms involved to one place and transcribe them into tidy readable/maintainable object-oriented c++ code. 

## Installation Notes
Running 'make' in the root directory should compile the source code a produce binary 'ccc_test' and shated library 'ccc.so'. The compilation was tested on Fedora 29 distribution.

## Main Components of the library

### Laguerre's Method
My code uses polynomial solving extensively. For that I use [Laguerre's Method](https://en.wikipedia.org/wiki/Laguerre%27s_method). My version of the algorithm is greatly inspired by book 'Numerical Recipes 3rd Edition: The Art of Scientific Computing' which is protected by Cambridge Press Copyright.
As I do not want to bother with the Copyright issues, I use compiled version of my algorithm which deviates from the published version. 

### CCC 
Class for producing Critical Curves and Caustics, which I use to characterise the microlenses. The class also includes a wrapper to be used within a python.

### main.cc
As this is mostly meant as a shared library to be later used within a Python application, main.cc serves for regression testing. Plan is to move to unit testing as the library grows larger. 

