%usefile c:\cprog\test\prg\nbinomial_poga_spline.prg

dataset d
d.infile using c:\cprog\test\testdata\nbinomial_poga_spline.raw

bayesreg b
b.outfile = c:\cprog\test\results\nbinomial_poga_spline
%b.regress resp = x(psplinerw2), iterations=12000 step=10 burnin=2000 family=nbinomial distopt=poga propopt=uniform predict using d

b.regress resp = x(psplinerw2), iterations=12000 step=10 burnin=2000 family=nbinomial distopt=nb propopt=uniform predict using d

b.getsample
b.autocor