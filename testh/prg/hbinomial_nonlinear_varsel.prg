% usefile c:\bayesx\trunk\testh\prg\hbinomial_nonlinear_varsel.prg

dataset d

d.infile using c:\bayesx\trunk\testh\testdata\hgaussian_nonlinear.raw
d.generate x3 = uniform()
d.generate x2 = uniform()

d.generate eta=0.45*sin(2*3.1415*x1)
d.generate pi = exp(eta)/(1+exp(eta))
d.generate ybin = binomial(1,pi)

d.outfile, header replace using c:\bayesx\trunk\testh\testdata\test_varsel_binomial.raw 

%d.infile using c:\bayesx\trunk\testh\testdata\test_varsel_binomial.raw

mcmcreg b

logopen using c:\bayesx\trunk\testh\results\binomial_nonlinear_varsel.log


%b.outfile = c:\bayesx\trunk\testh\results\binomial_nonlinear
%b.hregress ybin = const + x1 + x1(pspline,centermethod=nullspace) + x2 + x2(pspline,centermethod=nullspace), iterations=12000 step=10 burnin=2000 family=binomial_logit predict=light setseed=1234 using d

b.outfile = c:\bayesx\trunk\testh\results\binomial_nonlinear_varsel
<<<<<<< .mine
b.hregress ybin = const + x1 + x1(pspline,centermethod=nullspace,prior=ssvs) + x2 + x2(ssvs,centermethod=nullspace,prior=ssvs), iterations=12000 step=10 burnin=2000 family=binomial_logit predict=light setseed=1234 using d
=======
b.hregress ybin = const + x1(ssvs) + x1(pspline,centermethod=nullspace,prior=ssvs) + x2(ssvs) + x2(pspline,centermethod=nullspace,prior=ssvs), iterations=12000 step=10 burnin=2000 family=binomial_logit predict=light setseed=1234 using d
%b.hregress ybin = const + x1(ssvs) + x2(ssvs), iterations=12000 step=10 burnin=2000 family=binomial_logit predict=light setseed=1234 using d
%b.hregress ybin = const + x1 + x2, iterations=12000 step=10 burnin=2000 family=binomial_logit predict=light setseed=1234 using d
>>>>>>> .r1648

%b.hregress ybin = const + x2 + x2(pspline,centermethod=nullspace,prior=ssvs), iterations=12000 step=10 burnin=2000 family=binomial_logit predict=light setseed=1234 using d

%b.hregress ybin = const + x1(pspline,lambda=100,prior=ssvs)+ x2(pspline,lambda=100,prior=ssvs)+ x3(pspline,lambda=100,prior=ssvs), iterations=12000 step=10 burnin=2000 family=binomial_logit predict=light setseed=1234 using d
%b.hregress ybin = const + x1(pspline,lambda=100), iterations=12000 step=10 burnin=2000 family=binomial_logit predict=light setseed=1234 using d

b.getsample

logclose