% usefile c:\bayesx\testh\prg\gaussianh_nonp1_alt.prg

dataset d1
d1.infile using c:\bayesx\testh\testdata\hgaussian_nonp1_re.raw

dataset d2
d2.infile using c:\bayesx\testh\testdata\hgaussian_nonp1.raw


bayesreg b


b.outfile = c:\bayesx\testh\results\hg1_alt
b.regress y1 = x1(psplinerw2) + x2(psplinerw2)+ x3(psplinerw2) + x4(psplinerw2) + id(random) ,  iterations=12000 step=10 burnin=2000 family=gaussian using d2

% b.regress y1 = x3(psplinerw2) ,  iterations=12000 step=10 burnin=2000 family=gaussian using d2



