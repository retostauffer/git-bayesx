%usefile c:\cprog\test\prg\nbinomial_poig_spline.prg

dataset d
d.infile using c:\cprog\test\testdata\nbinomial_poig_spline.raw

bayesreg b
b.outfile = c:\cprog\test\results\nbinomial_poig_spline
b.regress resp = x(psplinerw2), iterations=12000 step=10 burnin=2000 family=nbinomial distopt=poig propopt=uniform predict using d

b.getsample
b.autocor