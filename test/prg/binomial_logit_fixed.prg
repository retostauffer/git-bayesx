% usefile c:\bayesx\test\prg\binomial_logit_fixed.prg

dataset d
d.infile using c:\bayesx\test\testdata\binomial_logit_fixed.raw

bayesreg b
b.outfile = c:\cprog\bayesx\results\binomial_logit_fixed
b.regress y = x1+x2 ,  predict   iterations=12000 step=10 burnin=2000 family=binomial using d