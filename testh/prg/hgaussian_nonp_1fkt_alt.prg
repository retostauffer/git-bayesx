% usefile c:\bayesx\testh\prg\hgaussian_nonp_1fkt_alt.prg

dataset d
d.infile using c:\bayesx\testh\testdata\hgaussian_nonp_1fkt.raw

bayesreg b
% b.outfile = c:\bayesx\testh\results\hgaussian_nonp_1fkt_alt
% b.regress y = x1(psplinerw2,contourprob=2) , predict iterations=12000 step=10 burnin=2000 family=gaussian using d


b.outfile = c:\bayesx\testh\results\hgaussian_nonp_1fkt_alt
b.regress y = x1(psplinerw2,contourprob=4) , predict iterations=12000 step=10 burnin=2000 family=gaussian using d
