%usefile c:\cprog\test\prg\nbinomial_poln_random.prg

dataset d
d.infile using c:\cprog\test\testdata\nbinomial_poln_random.raw
d.generate unit=_n

bayesreg b
b.outfile = c:\cprog\test\results\nbinomial_poln_random
b.regress resp = ind(random, lambda=1000)+unit(random, lambda=1000), iterations=12000 step=10 burnin=2000 family=poisson predict using d


b.getsample
b.autocor