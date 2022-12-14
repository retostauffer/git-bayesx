# BayesX version 3.0.2 ###

* Last changes: 9.3.2019

## Introduction

This Readme file provides information on how to compile BayesX from source code.

If you do not have the source code avalailable already, it can be obtained in two 
different ways:

- The download website on www.bayesx.org offers a zip-archive with the current
  stable release of BayesX (sourcecode.zip).
- There is an anonymous (read-only) SVN access to the current development version
  of BayesX via https://svn.gwdg.de/svn/bayesx/. This is not considered a stable
  distribution and may comprise development versions that may fail for several
  reasons.

For compiling the source code, we provide support for two different versions:

- Customized make support via the cmake functionality (http://www.cmake.org/).
- A basic makefile that has to be adapted manually to your OS.

The advantage of cmake is that it provides automated support for adapting the
makefile to your local computing environment (automatic detection of some
required components) across a variety of platforms. The ordinary makefile makes
less requirement for your system but may require more manual adujstments.

Comments / additions / corrections to the information in this file should be
submitted to bayesx@gwdg.de.


## Build BayesX using cmake

The basic idea of cmake is to generate a makefile for BayesX adapted to your
local working environment. In its somplest form, this should work as follows:

```
cd <location of the source code>
cmake .
```

This will generate the Makefile so that you can call

```
make
```

and this should finally yield the compiled BayesX versions (more details on the
two versions are given below). Unfortunately, this basic version will most
probably only work on platforms with a good preset of tools (e.g. most linux
distributions) but in most cases some additional arguments will be required (see
below).

The cmake specifications are made in CMakeLists.txt which comprise a specification 
of the sources included and additional information on required libraries etc. Usually, 
you should not have to make any changes to these files.

Calling cmake will create a make file that builds the BayesX executable, such that

```
bayesx
```

starts the software. Note that BayesX has to be executable. Under linux
systems, the call sometimes has to be changed to

```
./bayesx
```

### Requirements:

To be able to use cmake to generate a customized Makefile, you will need the
following components:

- the cmake (>= version 2.8.10) facility (available from www.cmake.org).
- C and C++ compiler such as gcc and g++ (available for windows via www.mingw.org).
- when working in a unix-type environment (linux, Mac OS), you will need the
  readline library. This is not required for windows platforms.

The directories containing cmake as well as the C and C++ compiler have to be in
your search path when executing the cmake command. For windows systems,
you can temporarily add a directory to the path by specifying for example

```
PATH c:\rtools\MinGW64\bin;c:\rtools\mingw\bin;%PATH%
```

which adds "c:\rtools\MinGW64\bin" and "c:\rtools\mingw\bin" (the locations of 32bit
and 64bit MinGW executables) to the current path. The command

```
echo %PATH%
```

displays the current path.

Note that paths including the directory of the r-tools binaries (such as c:\rtools\bin)
will conflict with cmake and should therefore be excluded from the path variable.

### Providing additional options for cmake

In most cases (at least under windows OS), additional options will have to be
parsed to cmake. Useful options from our experience are the following:

To explicitly specify the make facility that will be used for compiling BayesX,
set the `-G <GENERATOR>` option. For example, to use mingw32-make, you will have to
specify `-G "MinGW Makefiles"` as an option, i.e. the call will be

```
cmake -G "MinGW Makefiles" .
```

Calling cmake without options will give you a list of supported make facilities.

The C and C++ compiler that will be used by make can be specified by
`-D CMAKE_C_COMPILER=<C COMPILER> abd -D CMAKE_CXX_COMPILER=<C++ Compiler>`
For example to specify 64bit MinGW compilers, the cmake call may look as follows:

```
cmake -D CMAKE_C_COMPILER=x86_64-w64-mingw32-gcc -D CMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++ .
```

Of course, all options can be combined.

### Some known problems / trouble shooting:

- cmake conflicts with r-tools: If you have installed r-tools and included r-tools
binaries in your search path, this will conflict with cmake. You will have to delete r-tools
from the path (see "Requirements" above).

- If cmake failes, it sometimes helps to delete the cmake cache file CMakeCache.txt

- When checking out BayesX via SVN, you may conflict with an older makefile already contained
in the SVN directory. This old file is lower case (makefile) while the makefile generated
by cmake is upper case (Makefile). If your system is case specific, the call to make
may use either of the two. You may either delete the lower case makefile or specify explicitly
that Makefile shall be used by

```
make -f Makefile
```

### Some examples of working cmake specifications:

#### 64bit windows machine:

We assume that we start with a "clean" path that does not contain cmake,
the compilers, the make facility and no r-tools. Then we proceed as follows:

Add cmake to the path:
```
PATH C:\Program Files (x86)\cmake 2.8\bin;%PATH%
```

Add mingw (obtained as a part of the CodeBlocks compiler suite) to the path to
make mingw32-make available:
```
PATH c:\Program Files (x86)\CodeBlocks\MinGW\bin;%PATH%
```

Call cmake with explicit compiler specification for 64bit comilation and requesting
MinGW makefiles:
```
cmake -G "MinGW Makefiles" -D CMAKE_C_COMPILER=x86_64-w64-mingw32-gcc -D CMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++ .
```

The compiler specification will be stored in the cmake cache. Therefore, once you
have successfully called cmake once, the call
```
cmake -G "MinGW Makefiles" .
```
will be sufficient since the required compiler info is extracted from the cache.

Now call make
```
mingw32-make
```

If you are working in an environment with multiple kernels, the no. of kernels
to be used by make can be specified by the -j option. Four example, four kernels
will be used with
```
mingw32-make -j4
```

BayesX can now be invoked by
```
BayesX
```


#### 32bit windows machine

We assume that we start with a "clean" path that does not contain cmake,
the compilers, the make facility and no r-tools. Then we proceed as follows:

Add cmake to the path:
```
PATH C:\Programme\cmake 2.8\bin;%PATH%
```

Add mingw (obtained as a part of the CodeBlocks compiler suite) to the path to
make gcc/g++ available:
```
PATH c:\Programme\CodeBlocks\MinGW\bin;%PATH%
```

Call cmake with explicit compiler specification and requesting MinGW makefiles:
```
cmake -G "MinGW Makefiles" -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ .
```

The compiler specification will be stored in the cmake cache. Therefore, once you
have successfully called cmake once, the call
```
cmake -G "MinGW Makefiles" .
```
will be sufficient since the required compiler info is extracted from the cache.

Add r-tools to make mingw32-make available:
```
PATH c:\rtools\bin;C:\Rtools\MinGW\bin;C:\Rtools\MinGW64\bin;%PATH%
```
and call make (with no. of kernels defined by -j):
```
mingw32-make -j2
```

BayesX can now be invoked by
```
BayesX
```


## Build BayesX using makefile.orig

A basic makefile for BayesX is provided in makefile.orig. If you want to use this
makefile without relying on the cmake facility, either rename makefile.orig to makefile
and invoke

```
make BayesX
```

or call

```
make -f makefile.orig BayesX
```

This should compile an executable BayesX in the current working directory that
can be started by

```
BayesX
```

Note that BayesX has to be executable. Under linux systems, the call sometimes
has to be changed to

```
./BayesX
```



