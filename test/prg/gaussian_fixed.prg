% usefile c:\bayesx\trunk\test\prg\gaussian_fixed.prg

dataset d
d.infile using c:\bayesx\trunk\test\testdata\gaussian_fixed.raw

bayesreg b
b.outfile = c:\bayesx\trunk\test\results\gaussian_fixed
b.regress y = x1+x2, predict iterations=12000 burnin=2000 step=10 family=gaussian using d

%mcmcreg m
%m.outfile = c:\bayesx\trunk\test\results\gaussian_fixed_super
%m.hregress y = const + x1 + x2, predict=light iterations=12000 burnin=2000 step=10 family=gaussian using d