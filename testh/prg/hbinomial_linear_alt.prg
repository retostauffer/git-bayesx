% usefile c:\bayesx\testh\prg\hbinomial_linear_alt.prg

dataset d
d.infile using c:\bayesx\testh\testdata\hbinomial_linear.raw

bayesreg b
b.outfile = c:\bayesx\testh\results\hbinomial_linear_alt
b.regress y = x1+x2 , predict iterations=12000 step=10 burnin=2000 family=binomial using d

