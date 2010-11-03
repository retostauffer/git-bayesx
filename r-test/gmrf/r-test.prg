% usefile c:/bayesx/r-test/gmrf/r-test.prg

dataset d
d.infile using c:/bayesx/r-test/gmrf/data.raw

map m1 
m1.infile using c:/bayesx/r-test/gmrf/kreisesim.bnd

map m2
m2.infile, graph using c:/bayesx/r-test/gmrf/kreisesim.gra

remlreg r
bayesreg b
stepwisereg s

% remlreg

r.outfile = c:/bayesx/r-test/gmrf/output/reml_gauss_gam1
r.regress norm1 = region(spatial, map=m1), family=gaussian using d

r.outfile = c:/bayesx/r-test/gmrf/output/reml_gauss_gam2
r.regress norm1 = region(spatial, map=m2), family=gaussian using d

r.outfile = c:/bayesx/r-test/gmrf/output/reml_gauss_gam3
r.regress norm1 = region(geospline, map=m1, nrknots=12), family=gaussian using d

r.outfile = c:/bayesx/r-test/gmrf/output/reml_gauss_gam3
r.regress norm1 = region(geokriging, map=m1, nrknots=100), family=gaussian using d


r.outfile = c:/bayesx/r-test/gmrf/output/reml_gauss_vcm1
r.regress norm1 = region(spatial, map=m1) + z1*region(spatial, map=m1), family=gaussian using d

r.outfile = c:/bayesx/r-test/gmrf/output/reml_gauss_vcm2
r.regress norm1 = region(spatial, map=m1) + z1 + z1*region(spatial, map=m1, center), family=gaussian using d

r.outfile = c:/bayesx/r-test/gmrf/output/reml_gauss_vcm3
r.regress norm1 = region(spatial, map=m1) + z1*region(geospline, map=m1, nrknots=12), family=gaussian using d

r.outfile = c:/bayesx/r-test/gmrf/output/reml_gauss_vcm4
%r.regress norm1 = region(geokriging, map=m1) + z1*region(geokriging, map=m1, nrknots=100), family=gaussian using d


r.outfile = c:/bayesx/r-test/gmrf/output/reml_bernoulli_logit_gam1
r.regress bernlogit1 = region(spatial, map=m1), family=binomial using d

r.outfile = c:/bayesx/r-test/gmrf/output/reml_bernoulli_logit_vcm1
r.regress bernlogit2 = region(spatial, map=m1) + z1*region(spatial, map=m1), family=binomial using d


r.outfile = c:/bayesx/r-test/gmrf/output/reml_bernoulli_probit_gam1
r.regress bernprobit1 = region(spatial, map=m1), family=binomialprobit using d

r.outfile = c:/bayesx/r-test/gmrf/output/reml_bernoulli_probit_vcm1
r.regress bernprobit2 = region(spatial, map=m1) + z1*region(spatial, map=m1), family=binomialprobit using d


r.outfile = c:/bayesx/r-test/gmrf/output/reml_bernoulli_cloglog_gam1
r.regress berncloglog1 = region(spatial, map=m1), family=binomialcomploglog using d

r.outfile = c:/bayesx/r-test/gmrf/output/reml_bernoulli_cloglog_vcm1
%r.regress berncloglog2 = region(spatial, map=m1) + z1*region(spatial, map=m1), family=binomialcomploglog using d


r.outfile = c:/bayesx/r-test/gmrf/output/reml_binomial_gam1
r.regress bin1 = region(spatial, map=m1) weight binw, family=binomial using d

r.outfile = c:/bayesx/r-test/gmrf/output/reml_binomial_vcm1
r.regress bin2 = region(spatial, map=m1) + z1*region(spatial, map=m1) weight binw, family=binomial using d


r.outfile = c:/bayesx/r-test/gmrf/output/reml_poisson_gam1
r.regress pois1 = region(spatial, map=m1), family=poisson using d

r.outfile = c:/bayesx/r-test/gmrf/output/reml_poisson_vcm1
%r.regress pois2 = region(spatial, map=m1) + z1*region(spatial, map=m1), family=poisson using d


% bayesreg



% stepwisereg

