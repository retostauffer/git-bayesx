% usefile c:\bayesx\testh\prg\hbinomial_linear.prg

dataset d
d.infile using c:\bayesx\testh\testdata\hbinomial_linear.raw

mcmcreg b
b.outfile = c:\bayesx\testh\results\hbinomial_linear
b.hregress y = x1+x2 , predict=full iterations=12000 step=10 burnin=2000 family=binomial_logit using d

