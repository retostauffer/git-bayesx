% usefile c:\bayesx\testh\prg\multinom.prg

dataset d
d.infile using c:\bayesx\testh\testdata\multinomprobit.raw

mcmcreg b
b.outfile = c:\bayesx\testh\results\c1
b.hregress y1 = const + x1 + x2 , family=multinom_probit equationtype=meanservant hlevel=1 iterations=12000 step=10 burnin=2000 using d

b.outfile = c:\bayesx\testh\results\c1
b.hregress y2 = const + x1 + x2 , family=multinom_probit  equationtype=mean hlevel=1  using d


% mcmcreg c
% c.outfile = c:\bayesx\testh\results\pr1
% c.hregress y1 = const + x1 + x2 , family=binomialprobit modeonly hlevel=1 iterations=12000 step=10 burnin=2000 using d