%usefile c:\cprog\test\prg\nbinomial_poln_fixed.prg

dataset d
d.infile using c:\cprog\test\testdata\nbinomial_poln_fixed.raw
d.generate unit=_n

bayesreg b
b.outfile = c:\cprog\test\results\nbinomial_poln_fixed
b.regress resp = x+xx+unit(random, lambda=1000), iterations=12000 step=10 burnin=2000 family=poisson predict using d


b.getsample
b.autocor