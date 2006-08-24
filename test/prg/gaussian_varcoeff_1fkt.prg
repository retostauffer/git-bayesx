% usefile c:\cprog\test\prg\gaussian_varcoeff_1fkt.prg

dataset d
d.infile using c:\cprog\test\testdata\gaussian_varcoeff_1fkt.raw

bayesreg b
b.outfile = c:\cprog\test\results\gaussian_varcoeff_1fkt
b.regress y = x2*x1(psplinerw2) , iterations=12000 step=10 burnin=2000 family=gaussian predict using d

%b.regress y = x2*x1(rw1,a=1.0,b=0.005,lambda=0.1), iterations=12000 step=10 burnin=2000 family=gaussian predict using d
%b.regress y = x2*x1(rw2,a=1.0,b=0.005,lambda=0.1), iterations=12000 step=10 burnin=2000 family=gaussian predict using d

%b.regress y = x2*x1(psplinerw1,a=1.0,b=0.005,lambda=0.1,nrknots=20,degree=3), iterations=12000 step=10 burnin=2000 family=gaussian predict using d
%b.regress y = x2*x1(psplinerw2,a=1.0,b=0.005,lambda=0.1,nrknots=20,degree=3), iterations=12000 step=10 burnin=2000 family=gaussian predict using d

%b.regress y = x2*x1(nspline,a=1.0,b=0.005,lambda=0.1,nrknots=20), iterations=12000 step=10 burnin=2000 family=gaussian predict using d
%b.regress y = x2*x1(sspline,a=1.0,b=0.005,lambda=0.1), iterations=12000 step=10 burnin=2000 family=gaussian predict using d
