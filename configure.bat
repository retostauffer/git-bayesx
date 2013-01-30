# 64 bit Windows (Thomas' office PC)
# add mingw to the path
# add cmake to the path
cmake -G "MinGW Makefiles" -D CMAKE_C_COMPILER=x86_64-w64-mingw32-gcc -D CMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++ .

# 32 bit Windows (Thomas' Notebook)
# add mingw to the path
# add cmake to the path
# if there has been no previous call that set the compiler:
cmake -G "MinGW Makefiles" -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ .
# otherwise:
cmake -G "MinGW Makefiles" .
