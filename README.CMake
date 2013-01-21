Building BayesX using CMake.
============================

CMake is a meta-build system that supports a large variety of build systems
and integrated development environments.
I recommend using CMake Version 2.8.7 or later.

The main file to control the 'products' to build are specified in the 
CMakeLists.txt files of the source tree.

At this stage we have two products, 
 - executable - defined in CMakeLists.txt and
 - shared linked library for the java version - defined in java/CMakeLists.txt.

When CMake is invoked, it generates build tool-specific files.
After that one starts the IDE or tool of choice and opens the project files.

Note that one never adds files to via the IDE tool but adds the files
to the CMakeLists.txt instead. The syntax is quite easy and one has
a lot of options for searching additional libs etc..

How to build BayesX with CMake:

  1. Step: Build the generator-specific files (e.g. GNU Makefile)
     
     # cmake -G <GENERATOR> <SOURCE-DIR>

      on most platforms with a good preset of tools (e.g. Linux, Mac OS X)
      the following works fine:
    
     # cmake .

  2. Step: Build the Project

     # make

     This will build bayesx executable and java/BayesXDll.


Building on Unix using 'make':
------------------------------

  # cmake .
  # make


Building on Windows:
-------------------- 

Using MinGW32 and MinGW32 Make:

  > cmake -G "MinGW Makefiles"
  > mingw32-make

Using CodeBlocks with MinGW32:

  > cmake -G "CodeBlocks - MinGW Makefiles"
  open CodeBlocks file "BayesX.cbc" Project file and hit build button.

Using 64-bit Compiler from RTools and Codeblocks:

  Put Rtools and Codeblocks into path.

  > cmake -G "MinGW Makefiles" -D CMAKE_C_COMPILER=x86_64-w64-mingw32-gcc -D CMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++ .
  > mingw32-make

  See also ./configure.bat batch file for this configuration.


Building on Mac OS X:
---------------------

Using GNU make on shell:

  # cmake .
  # make

Using Xcode:

  # cmake -G "Xcode"
  open BayesX.xcodeproj and hit build button.


Run BayesX Java GUI
===================

  1. cd java
  2. javac BayesX.java
  3. java BayesX


Adding new sources
==================

New file names should be added to CMake variables in CMakeLists.txt and/or 
java/CMakeLists.txt. You see them at the very beginning. Currently they
are splitted (for historical reason) into various categories.
See 'add_executable' and 'add_library' where they are ultimately assembled
together. 
 

Explicitly build only executable or only dll
============================================

make <target> 
  where target is bayesx or BayesXDll


Build faster on multi-processor
===============================

make -j#

  where # is the number of cores (e.g. make -j4)


Last Updated: 2013-01-21 <dadler@uni-goettingen.de>

