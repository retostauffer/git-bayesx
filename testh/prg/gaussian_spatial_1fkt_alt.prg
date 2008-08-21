% usefile c:\bayesx\testh\prg\gaussian_spatial_1fkt_alt.prg

dataset d
d.infile using c:\bayesx\testh\testdata\gaussian_spatial_1fkt.raw

map m
m.infile using c:\bayesx\testh\testdata\kreisesim.bnd
m.reorder

bayesreg b
b.outfile = c:\bayesx\testh\results\gaussian_spatial_1fkt_alt
b.regress y1 = region(spatial,map=m,lambda=10) , modeonly family=gaussian iterations=12000 step=10 burnin=2000 using d

%b.regress y1 = region(geospline,map=m,a=1.0,b=0.005,lambda=0.1,degree=3,nrknots=20), family=gaussian step=10 burnin=2000 iterations=12000 predict using d
