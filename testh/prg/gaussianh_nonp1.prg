% usefile c:\bayesx\testh\prg\gaussianh_nonp1.prg

dataset d1
d1.infile using c:\bayesx\testh\testdata\hgaussian_nonp1_re.raw

dataset d2
d2.infile using c:\bayesx\testh\testdata\hgaussian_nonp1.raw


mcmcreg b

b.outfile = c:\bayesx\testh\results\hg2_1
b.hregress id = x1(pspline,centermethod=nullspace) + x2(pspline,centermethod=nullspace), iterations=12000 step=10 burnin=2000 family=gaussian_re hlevel=2 using d1

b.outfile = c:\bayesx\testh\results\hg1_1
b.hregress y1 =  x3(pspline,centermethod=nullspace) + x4(pspline,centermethod=nullspace) + id(hrandom)  ,  iterations=12000 step=10 burnin=2000 family=gaussian hlevel=1 using d2




