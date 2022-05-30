# The Channel2D Libray  

Code for threading uniform, cubic b-spline curves into polygonsl channels.

## INTRODUCTION

The code is organized as follows:

* bin         - subdirectory where the executable channel2d-app will be written to     
* data        - subdirectory where the example input files and empty output subdirectories are 
* doc         - subdirectory where doxygen documentation files are
* include     - subdirectory where include files for the library channel2d will be written to
* lib         - subdirectory where the lib file for the library channel2d will be written to
* LICENSE.md  - license file
* README.md   - this file
* script      - subdirectory containing a shell script to execute the examples
* src         - subdirectory containing the source files of the library channel2d and application channe2d-app

A detailed documentation of the code can be found in the file refman.pdf inside directory doc.

In a nutshell, this code implements a planar version of the algorithm described at

Ashish Myles and Jorg Peters, Threading splines through 3D channels, Computer-Aided Design, 37(2), 2005, p. 139-148.

## INSTALLATION

Please, take a look at Chapter 3 of the document \c refman.pdf inside subdirectory \c doc.

## EXAMPLES

You can run the executable 'channel2d-app' on the data files in subdirectory  

  <path to directory channel2d>/data  

by executing the shell script 'run.sh' that you find in subdirectory  

  <path to directory channel2d>/script  

##  LAST UPDATE

May 29, 2022

## CONTACT

If you run  into trouble compiling or using the library, please email me at:

mfsiqueira@gmail.com

Have fun!

Marcelo Siqueira
