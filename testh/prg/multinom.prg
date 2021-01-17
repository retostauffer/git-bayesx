% usefile c:\bayesx\trunk\testh\prg\multinom.prg

dataset d
d.infile using c:\bayesx\trunk\testh\testdata\multinomprobit.raw

mcmcreg b
b.outfile = c:\bayesx\trunk\testh\results\c1
b.hregress y1 = const + x1 + x2 , family=multinom_probit equationtype=servant hlevel=1 iterations=12000 step=10 burnin=2000 using d

b.outfile = c:\bayesx\trunk\testh\results\c2
b.hregress y2 = const + x1 + x2 , family=multinom_probit equationtype=main hlevel=1  using d


% mcmcreg c
% c.outfile = c:\bayesx\trunk\testh\results\pr1
% c.hregress y1 = const + x1 + x2 , family=binomialprobit modeonly hlevel=1 iterations=12000 step=10 burnin=2000 using d