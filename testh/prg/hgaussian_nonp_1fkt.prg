% usefile c:\bayesx\testh\prg\hgaussian_nonp_1fkt.prg

dataset d
d.infile using c:\bayesx\testh\testdata\hgaussian_nonp_1fkt.raw

mcmcreg b
b.outfile = c:\bayesx\testh\results\hgaussian_nonp_1fkt
b.hregress y = x1(pspline,meaneffect) ,   iterations=12000 step=10 burnin=2000 family=gaussian using d


% b.outfile = c:\bayesx\testh\results\hgaussian_nonp_1fkt
% b.hregress y = x1(pspline,constraints=increasing) ,   iterations=12000 step=10 burnin=2000 family=gaussian using d
