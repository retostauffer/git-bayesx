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
cp -rf examples sourcecode

cd sourcecode\adaptiv
rm -rf CVS
rm *.~*
rm #*
rm *.obj

cd ..\examples
rm -rf CVS
rm *.~*
rm #*
rm *.obj

cd ..\alex
rm -rf CVS
rm *.~*
rm #*
rm *.obj

cd ..\andrea
rm -rf CVS
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
rm -rf CVS
rm *.~*
rm #*
rm *.obj

cd ..\graph
rm -rf CVS
rm *.~*
rm #*
rm *.obj

cd ..\leyre
rm -rf CVS
rm *.~*
rm #*
rm *.obj

cd ..\mcmc
rm -rf CVS
rm *.~*
rm #*
rm *.obj

cd ..\psplines
rm -rf CVS
rm *.~*
rm #*
rm *.obj

cd ..\samson
rm -rf CVS
rm *.~*
rm #*
rm *.obj

cd ..\structadd
rm -rf CVS
rm *.~*
rm #*
rm *.obj

cd ..
mkdir gnuobj

rm bayesxsource.zip
zip -ll -r bayesxsource.zip c:\bayesx\sourcecode


