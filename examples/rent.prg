dataset rent
rent.infile using c:\bayes\examples\rent94.raw

bayesreg b
b.outfile = c:\results\rent
b.regress R = const + F(rw2) + A(rw2) + L(random) , iterations=22000 burnin=2000 step=20 family=gaussian using rent
 
b.outfile = c:\results\rent_2
b.regress R = const + F(rw2) + A(rw2) + L(random) if A > 1890 , iterations=22000 burnin=2000 step=20 family=gaussian using rent