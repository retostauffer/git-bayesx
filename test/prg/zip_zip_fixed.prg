%usefile c:\cprog\test\prg\nbinomial_poga_fixed.prg

dataset d
d.infile using c:\cprog\test\testdata\nbinomial_poga_fixed.raw

bayesreg b
b.outfile = c:\cprog\test\results\nbinomial_poga_fixed
b.regress resp = x+xx, iterations=12000 step=10 burnin=2000 family=nbinomial distopt=poga propopt=uniform predict using d

%b.regress resp = x+xx, iterations=12000 step=10 burnin=2000 family=nbinomial distopt=nb propopt=uniform predict using d

b.getsample
b.autocor