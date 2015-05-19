% usefile c:\bayesx\trunk\testh\prg\hbinomial_nonlinear_varsel.prg

dataset d
d.infile using c:\bayesx\trunk\testh\testdata\hgaussian_nonlinear.raw
d.generate x3 = uniform()
d.generate x2 = uniform()

d.generate eta=sin(2*3.1415*x1)
d.generate pi = exp(eta)/(1+exp(eta))
d.generate ybin = binomial(1,pi)


%d.outfile, header replace using c:\tmp\test.raw 

mcmcreg b

logopen using c:\bayesx\trunk\testh\results\binomial_nonlinear_varsel.log

b.outfile = c:\bayesx\trunk\testh\results\binomial_nonlinear_varsel

%b.hregress ybin = const + x1(pspline,lambda=100,prior=ssvs)+ x2(pspline,lambda=100,prior=ssvs)+ x3(pspline,lambda=100,prior=ssvs), iterations=12000 step=10 burnin=2000 family=binomial_logit predict=light using d
b.hregress ybin = const + x1(pspline,lambda=100), iterations=12000 step=10 burnin=2000 family=binomial_logit predict=light using d

b.getsample

logclose