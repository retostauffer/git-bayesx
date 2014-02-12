% usefile c:\lab\bayes_commit\trunk\testh\prg\hgaussian_linear.prg

dataset d
d.infile using c:\tmp\d.raw
d.generate c = 1



mcmcreg b

b.outfile = c:\lab\bayes_commit\trunk\testh\results\hgaussian_linear
b.hregress y = c ,  iterations=12000 step=10 burnin=2000 family=hetgaussian equationtype = variance using d

b.outfile = c:\lab\bayes_commit\trunk\testh\results\hgaussian_linear
b.hregress y = c(ridge) ,  iterations=12000 step=10 burnin=2000 family=hetgaussian equationtype = mean using d

