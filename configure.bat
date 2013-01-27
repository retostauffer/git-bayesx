# 64 bit Windows (Thomas PC im Buero)
cmake -G "MinGW Makefiles" -D CMAKE_C_COMPILER=x86_64-w64-mingw32-gcc -D CMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++ .

# 32 bit Windows (Thomas Notebook)
cmake -G "MinGW Makefiles" -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ .
