%usefile c:\cprog\test\prg\nbinomial_poig_fixed.prg

dataset d
d.infile using c:\cprog\test\testdata\nbinomial_poig_fixed.raw

bayesreg b
b.outfile = c:\cprog\test\results\nbinomial_poig_fixed
b.regress resp = x+xx, iterations=12000 step=10 burnin=2000 family=nbinomial distopt=poig propopt=uniform predict using d

b.getsample
b.autocor