%usefile c:\cprog\test\prg\nbinomial_poga_random.prg

dataset d
d.infile using c:\cprog\test\testdata\nbinomial_poga_random.raw

bayesreg b
b.outfile = c:\cprog\test\results\nbinomial_poga_random
b.regress resp = ind(random), iterations=12000 step=10 burnin=2000 family=nbinomial distopt=poga propopt=uniform predict using d

%b.regress resp = ind(random), iterations=12000 step=10 burnin=2000 family=nbinomial distopt=nb propopt=uniform predict using d

b.getsample
b.autocor