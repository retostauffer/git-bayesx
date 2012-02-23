Building BayesX using CMake.
============================

CMake is a meta-build system that supports a large variety of build systems,
integrated development environments.
I recommend using CMake Version 2.8.7.

The main file to control the 'products' to build are specified in the 
CMakeLists.txt file at the top-level of the source tree.
Currently we have two products, the executable and a shared library.
 
When CMake is invoked, it generates build tool-specific files.
After that one starts the IDE or tool of choice and opens the project files.

Note, that one never adds files to via the IDE tool but adds the files
to the CMakeLists.txt instead. The syntax is quite easy and one has
a lot of options for searching additional libs etc..

Enjoy,
- Daniel



Building on Unix:
-----------------

Using GNU Make:
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


Building on Mac OS X:
---------------------

Using GNU make on shell:

  # cmake .
  # make

Using Xcode:

  # cmake -G "Xcode"
  open BayesX.xcodeproj and hit build button.



