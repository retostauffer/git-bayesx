% usefile c:\temp\hgaussian_nonlinear_varsel.prg

dataset d
d.infile using c:\temp\hgaussian_nonlinear.raw
d.generate x3 = uniform()
d.generate x2 = uniform()

d.generate y2 = 1*(uniform()<0.5)

d.descriptive y y2 x1 x2 x3

d.outfile, header replace using c:\temp\test.raw

mcmcreg b

logopen using c:\temp\hgaussian_nonlinear_varsel.log

b.outfile = c:\temp\hgaussian_nonlinear_varsel

%b.hregress y2 = const + x1(pspline,prior=ssvs) + x2(pspline,prior=ssvs), iterations=12000 step=10 burnin=2000 family=binomial_logit using d
%b.hregress y2 = const + x2(pspline), iterations=12000 step=10 burnin=2000 family=binomial_logit using d
%b.getsample

logclose