% usefile c:\cprog\test\prg\binomial_logit_fixed.prg

dataset d
d.infile using c:\cprog\test\testdata\binomial_logit_fixed.raw

bayesreg b
b.outfile = c:\cprog\test\results\binomial_logit_fixed
b.regress y = x1+x2 ,  predict  modeonly iterations=12000 step=10 burnin=2000 family=binomial using d