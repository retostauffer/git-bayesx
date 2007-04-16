% usefile c:\bayesx\test\prg\gaussian_mult_work.prg

dataset d
d.infile using c:\bayesx\test\testdata\gaussian_mult.raw

bayesreg b
b.outfile = c:\bayesx\test\results\gaussian_mult
b.regress y = cl*x(random_rw2) , modeonly iterations=12000 step=10 burnin=2000 family=gaussian predict using d

