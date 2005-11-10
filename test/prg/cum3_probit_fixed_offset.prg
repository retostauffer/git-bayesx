% usefile c:\bayesx\test\prg\cum3_probit_fixed_work.prg

dataset d
d.infile using c:\bayesx\test\testdata\cum3_probit_fixed_sort.raw


bayesreg b
b.outfile = c:\bayesx\test\results\cum3_probit_fixed
b.regress y = x1+x2 , predict iterations=12000 step=10 burnin=2000 family=cumprobit nosort using d

dataset s
s.infile using c:\bayesx\test\results\cum3_probit_fixed_predictmean.raw
s.generate o = linpred*sqrt(0.69)

b.regress y = o(offset)+x1+x2 , predict iterations=12000 step=10 burnin=2000 family=cumprobit nosort using s
