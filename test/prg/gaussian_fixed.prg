% usefile c:\bayesx\test\prg\gaussian_fixed.prg

dataset d
d.infile using c:\bayesx\test\testdata\gaussian_fixed.raw

bayesreg b
b.outfile = c:\bayesx\test\results\gaussian_fixed
b.regress y = x1+x2 , predict iterations=12000 burnin=2000 step=10 family=gaussian using d