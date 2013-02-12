cp -f CMakeLists.txt sourcecode\CMakeLists.txt
cp -f README.BayesX sourcecode\README.BayesX
cp -f makefile.orig sourcecode\makefile.orig
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
cp -rf share sourcecode

cp -f internet\bayesx\tutorials\zambia.bnd sourcecode\examples\zambia.bnd
cp -f internet\bayesx\tutorials\zambia.raw sourcecode\examples\zambia.raw
cp -f internet\bayesx\tutorials\regression_mcmc\mcmctutorial.prg sourcecode\examples\mcmctutorial.prg
cp -f internet\bayesx\tutorials\regression_reml\remltutorial.prg sourcecode\examples\remltutorial.prg
cp -f internet\bayesx\tutorials\regression_step\steptutorial.prg sourcecode\examples\steptutorial.prg

cd sourcecode\adaptiv
rm -rf CVS
rm -rf .svn
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\examples
rm -rf CVS
rm -rf .svn
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\alex
rm -rf CVS
rm -rf .svn
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\andrea
rm -rf CVS
rm -rf .svn
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\bib
rm -rf .svn
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm -f *.dfm
rm *.o
rm *.d

cd ..\dag
rm -rf .svn
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\graph
rm -rf .svn
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\leyre
rm -rf CVS
rm -rf .svn
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\mcmc
rm -rf CVS
rm -rf .svn
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\psplines
rm -rf CVS
rm -rf .svn
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\samson
rm -rf CVS
rm -rf .svn
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\structadd
rm -rf .svn
rm -rf CVS
rm *.~*
rm #*
rm *.obj
rm *.o
rm *.d

cd ..\share
rm -rf .svn
rm -rf CVS

cd cmake
rm -rf .svn
rm -rf CVS
cd ..

cd ..
mkdir java

cd ..
cp -rf java\doc sourcecode\java\doc
cp -rf java\umontreal sourcecode\java\umontreal

cd java
cp *.java ..\sourcecode\java
cp *.cpp ..\sourcecode\java
cp *.h ..\sourcecode\java
cp *.sh ..\sourcecode\java
cp *.txt ..\sourcecode\java
cp *.gif ..\sourcecode\java
cp *.jpg ..\sourcecode\java

cd ..\sourcecode
mkdir gnuobj

rm bayesxsource.zip
rm *.o
rm *.d
rm *.exe
zip -ll -r bayesxsource.zip .
zip -r bayesxsource.zip java\umontreal java\doc java\Bayesicon.gif

cd ..

