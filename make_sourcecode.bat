c:
cd c:\bayesx

cp -f makefile sourcecode\makefile
cp -f values.h sourcecode\values.h
cp -f main.cpp sourcecode\main.cpp

cp -rf adaptiv sourcecode
cp -rf alex sourcecode
cp -rf andrea sourcecode
cp -rf bib sourcecode
cp -rf dag sourcecode
cp -rf graph sourcecode
cp -rf leyre sourcecode
cp -rf mcmc sourcecode
cp -rf psplines sourcecode
cp -rf samson sourcecode
cp -rf structadd sourcecode

cd sourcecode\adaptiv
rm -r CVS
rm *.~*
rm #*
rm *.obj

cd ..\alex
rm -r CVS
rm *.~*
rm #*
rm *.obj

cd ..\andrea
rm -r CVS
rm *.~*
rm #*
rm *.obj

cd ..\bib
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm -f *.dfm

cd ..\dag
rm -r CVS
rm *.~*
rm #*
rm *.obj

cd ..\graph
rm -r CVS
rm *.~*
rm #*
rm *.obj

cd ..\leyre
rm -r CVS
rm *.~*
rm #*
rm *.obj

cd ..\mcmc
rm -r CVS
rm *.~*
rm #*
rm *.obj

cd ..\psplines
rm -r CVS
rm *.~*
rm #*
rm *.obj

cd ..\samson
rm -r CVS
rm *.~*
rm #*
rm *.obj

cd ..\structadd
rm -r CVS
rm *.~*
rm #*
rm *.obj

cd ..
mkdir gnuobj

rm bayesxsource.zip
zip -ll -r bayesxsource.zip c:\bayesx\sourcecode
