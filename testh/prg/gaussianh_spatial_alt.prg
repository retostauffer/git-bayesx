% usefile c:\bayesx\testh\prg\gaussianh_spatial_alt.prg


dataset d
d.infile using c:\bayesx\testh\testdata\gaussianh_spatial.raw

map m
m.infile using c:\bayesx\testh\testdata\kreisesim.bnd
m.reorder

delimiter = ;

bayesreg b;
b.outfile = c:\bayesx\testh\results\gaussianh_spatial_alt;
b.regress y1 =  region(random,lambda=10)+
region(spatial,map=m,lambda=10) +x1(psplinerw2)+x2(psplinerw2), iterations=12000 
step=10 burnin=2000 family=gaussian using d;


% b.regress y1 = x1(psplinerw2,lambda=10) + region(random,lambda=10)+
% region(spatial,map=m,lambda=10) , iterations=12000 
% step=10 burnin=2000 family=gaussian modeonly constlambda using d;


