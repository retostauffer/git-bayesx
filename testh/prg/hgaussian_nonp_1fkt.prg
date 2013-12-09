% usefile c:\bayesx\trunk\testh\prg\hgaussian_nonp_1fkt.prg

dataset d
d.infile using c:\bayesx\trunk\testh\testdata\hgaussian_nonp_1fkt.raw
d.generate x2 = normal()

mcmcreg b
b.outfile = c:\bayesx\trunk\testh\results\hgaussian_nonp_1fkt
b.hregress y = const + x1(pspline) ,   iterations=12000 step=10 burnin=2000 family=gaussian using d


% b.outfile = c:\bayesx\trunk\testh\results\hgaussian_nonp_1fkt
% b.hregress y = x1(pspline,constraints=increasing) ,   iterations=12000 step=10 burnin=2000 family=gaussian using d
