% usefile c:\bayesx\trunk\testh\prg\hgaussian_nonlinear_varsel.prg

dataset d
d.infile using c:\bayesx\trunk\testh\testdata\hgaussian_nonlinear.raw
d.generate x3 = uniform()
d.generate x2 = uniform()

%d.outfile, header replace using c:\tmp\test.raw 

mcmcreg b

logopen using c:\bayesx\trunk\testh\results\hgaussian_nonlinear_varsel.log

b.outfile = c:\bayesx\trunk\testh\results\hgaussian_nonlinear_varsel

b.hregress y = const + x1(pspline,lambda=100,prior=ssvs)+ x2(pspline,lambda=100,prior=ssvs)+ x3(pspline,lambda=100,prior=ssvs) , iterations=12000 step=10 burnin=2000 family=gaussian  using d

b.getsample

logclose