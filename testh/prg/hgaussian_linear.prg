% usefile c:\bayesx\testh\prg\hgaussian_linear.prg

dataset d
d.infile using c:\bayesx\testh\testdata\hgaussian_linear.raw

mcmcreg b
b.outfile = c:\bayesx\testh\results\hgaussian_linear
b.hregress y = x1+x2 ,  predict=full iterations=12000 step=10 burnin=2000 family=gaussian using d

