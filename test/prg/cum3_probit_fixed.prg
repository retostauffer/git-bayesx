% usefile c:\bayesx\test\prg\cum3_probit_fixed.prg

dataset d
d.infile using c:\bayesx\test\testdata\cum3_probit_fixed.raw

bayesreg b
b.outfile = c:\bayesx\test\results\cum3_probit_fixed
b.regress y = x1+x2 , predict iterations=12000 step=10 burnin=2000 family=cumprobit  using d
