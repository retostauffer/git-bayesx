% usefile c:\cprog\test\prg\nbinomial_poig_random.prg

dataset d
d.infile using c:\cprog\test\testdata\nbinomial_poig_random.raw

bayesreg b
b.outfile = c:\cprog\test\results\nbinomial_poig_random
b.regress resp = ind(random), iterations=12000 step=10 burnin=2000 family=nbinomial distopt=poig propopt=uniform predict using d

b.getsample
b.autocor