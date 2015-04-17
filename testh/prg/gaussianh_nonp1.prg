% usefile c:\bayesx\trunk\testh\prg\gaussianh_nonp1.prg

dataset d1
d1.infile using c:\bayesx\trunk\testh\testdata\hgaussian_nonp1_re.raw

dataset d2
d2.infile using c:\bayesx\trunk\testh\testdata\hgaussian_nonp1.raw

logopen c:\bayesx\trunk\testh\results\test.log

mcmcreg b

b.outfile = c:\bayesx\trunk\testh\results\hg2_1
b.hregress id = x1(pspline) + x2(pspline), iterations=12000 step=10 burnin=2000 family=gaussian_re hlevel=2 using d1

b.outfile = c:\bayesx\trunk\testh\results\hg1_1
b.hregress y1 =  x3(pspline) + x4(pspline) + id(hrandom)  ,  iterations=52000 step=50 burnin=2000 family=gaussian hlevel=1 using d2




