% usefile c:\bayesx\testh\prg\gaussianh_spatial_work.prg

dataset d
d.infile using c:\bayesx\testh\testdata\gaussianh_spatial_cl.raw
d.generate region2 = region

dataset d2
d2.infile using c:\bayesx\testh\testdata\gaussianh_spatial.raw

map m
m.infile using c:\bayesx\testh\testdata\kreisesim.bnd
m.reorder


mcmcreg b
b.outfile = c:\bayesx\testh\results\gaussianh_spatial
b.hregress region = x1(pspline,centermethod=nullspace)+region2(spatial,map=m) , hlevel=2 iterations=12000 step=10 burnin=2000 family=gaussian_re using d
b.hregress y1 =  x2(pspline,centermethod=nullspace)+region(hrandom)  , hlevel=1 iterations=12000 step=10 burnin=2000 family=gaussian using d2


% bayesreg b
% b.outfile = c:\bayesx\testh\results\gaussianh_spatial
% b.hregress region = region2(spatial,map=m,lambda=10) , constlambda modeonly iterations=12000 step=10 burnin=2000 family=gaussian_re using d
% b.hregress y1 = x1(psplinerw2,lambda=10) + region(hrandom,lambda=10)  , constlambda  modeonly iterations=12000 step=10 burnin=2000 family=gaussian using d2


