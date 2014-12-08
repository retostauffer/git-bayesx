% usefile c:\bayesx\trunk\testh\prg\hgaussian_linear.prg

dataset d
d.infile using c:\bayesx\trunk\testh\testdata\hgaussian_linear.raw
d.generate x3 = uniform()


mcmcreg b

logopen using c:\bayesx\trunk\testh\results\test.log

b.outfile = c:\bayesx\trunk\testh\results\hgaussian_linear
%b.hregress y = const + x3(ssvs,a=5,b=25), iterations=12000 step=10 burnin=2000 family=gaussian  using d
b.hregress y = const + x1(ssvs,a=5,b=25) + x2(ssvs,a=5,b=25) + x3(ssvs,a=5,b=25) , iterations=12000 step=10 burnin=2000 family=gaussian  using d


%  

