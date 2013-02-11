######################################
# 64 bit Windows (Thomas' office PC) #
######################################

# add r-tools mingw to the path to enable access to the compiler 
c:\env\mingw-rtools.bat
# add cmake to the path 
c:\env\cmake.bat
# add codeblocks mingw to the path to enable access to mingw32-make 
c:\env\mingw-cb.bat

# call cmake with explicit compiler specification for 64bit comilation
# if there has been no previous call that set the compiler:
cmake -G "MinGW Makefiles" -D CMAKE_C_COMPILER=x86_64-w64-mingw32-gcc -D CMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++ .
#otherwise
cmake -G "MinGW Makefiles" .

# call make (with no. of kernels defined by -j)
mingw32-make -j4

# add java to the path
c:\env\java.bat
# compile java version
cd java
javac BayesX.java
java BayesX


#####################################
# 32 bit Windows (Thomas' Notebook) #
#####################################

# add mingw to the path 
c:\env\mingw-cb.bat
# add cmake to the path 
c:\env\cmake.bat

# if there has been no previous call that set the compiler:
cmake -G "MinGW Makefiles" -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ .
# otherwise:
cmake -G "MinGW Makefiles" .

# make mingw32-make available
c:\env\rtools.bat
# call make (with no. of kernels defined by -j)
mingw32-make -j2

# make java available
c:\env\java.bat
# compile java version
cd java
javac BayesX.java
java BayesX
