c:
cd c:\bayesx

cp -f makefile sourcecode\makefile
cp -f values.h sourcecode\values.h
cp -f export_type.h sourcecode\export_type.h
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
rm *.o
rm *.d

cd ..\examples
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\alex
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\andrea
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\bib
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm -f *.dfm
rm *.o
rm *.d

cd ..\dag
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\graph
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\leyre
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\mcmc
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\psplines
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\samson
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\structadd
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..
mkdir gnuobj

rm bayesxsource.zip
rm *.o
rm *.d
rm *.exe
zip -ll -r bayesxsource.zip c:\bayesx\sourcecode


