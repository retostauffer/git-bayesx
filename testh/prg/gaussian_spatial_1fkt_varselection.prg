% usefile c:\bayesx\trunk\testh\prg\gaussian_spatial_1fkt_varselection.prg

dataset d
d.infile using c:\bayesx\trunk\testh\testdata\gaussian_spatial_1fkt.raw
d.generate y2 = normal()

map m
m.infile using c:\bayesx\trunk\testh\testdata\kreisesim.bnd


mcmcreg b
b.outfile = c:\bayesx\trunk\testh\results\gaussian_spatial_1fkt_work
b.hregress y1 = const + region(spatial,map=m,centermethod=nullspace,prior=ssvs), family=gaussian iterations=12000 step=10 burnin=2000 using d
b.hregress y2 = const + region(spatial,map=m,centermethod=nullspace,prior=ssvs), family=gaussian iterations=12000 step=10 burnin=2000 using d

