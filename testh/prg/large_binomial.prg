% usefile c:/bayesx/trunk/testh/prg/large_binomial.prg

dataset d
d.infile using c:/bayesx/trunk/testh/testdata/large_binomial.raw

mcmcreg m
m.outfile = c:/bayesx/trunk/testh/results/large_binomial
m.hregress y = x1(pspline, nocenter) + x2(pspline) + x3(pspline) + x4(pspline) + x5(pspline), family=binomial_logit iterations=12000 burnin=2000 step=10 setseed=5678 highspeedon using d
m.getsample

% drop d m

 
