% usefile c:\bayesx\test\prg\gaussian_spatial_vcm_1fkt_center.prg

dataset d
d.infile using c:\bayesx\test\testdata\gaussian_spatial_vcm_center.raw

map m
m.infile using c:\bayesx\test\testdata\kreisesim.bnd
m.reorder

bayesreg b
b.outfile = c:\bayesx\test\results\gaussian_spatial_vcm_1fkt
b.regress y1 = x2+x2*region(geospline,map=m,a=0.001,b=0.001,center,lambda=1000) , modeonly predict family=gaussian iterations=12000 step=10 burnin=2000 using d

stepwisereg s
 s.outfile = c:\bayesx\test\results\gaussian_varcoeff_1fkt_s
 s.regress y1 = x2(factor,coding=effect,reference=-1) + x2*region(geospline,map=m,center,sp,spmin=1000,spmax=2000,number=2,spstart=1000) , startmodel=userdefined steps=0 family=gaussian using d
% s.outfile = c:\bayesx\test\results\gaussian_varcoeff_1fkt_s
% s.regress y1 = x2(factor,coding=effect,reference=-1) + x2*region(spatial,map=m,center) ,  family=gaussian using d
