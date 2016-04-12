% usefile c:\bayesx\trunk\testh\prg\hbinomial_linear.prg

dataset d
d.infile using c:\bayesx\trunk\testh\testdata\hbinomial_linear.raw

logopen using c:\bayesx\trunk\results\hbinomial_linear.log

mcmcreg b
b.outfile = c:\bayesx\trunk\testh\results\hbinomial_linear
b.hregress y = const+x1+x2(pspline) , predict=light  iterations=12000 step=10 burnin=2000 family=binomial_logit using d

b.getsample