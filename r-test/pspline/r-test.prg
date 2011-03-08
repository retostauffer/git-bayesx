% usefile c:/bayesx/r-test/pspline/r-test.prg

dataset d
d.infile using c:/bayesx/r-test/pspline/data.raw

remlreg r
bayesreg b
stepwisereg s

% remlreg

r.outfile = c:/bayesx/r-test/pspline/output/reml_gauss_gam1
r.regress norm1 = x1(psplinerw2) + x2(psplinerw2), family=gaussian using d

r.outfile = c:/bayesx/r-test/pspline/output/reml_gauss_gam2
r.regress norm1 = x1(psplinerw2, degree=2, nrknots=12) + x2(psplinerw1, degree=1, nrknots=30), family=gaussian using d

r.outfile = c:/bayesx/r-test/pspline/output/reml_gauss_gam3
r.regress norm1 = x1(psplinerw1, degree=0, nrknots=12) + x2(psplinerw2, degree=1, nrknots=30), family=gaussian using d

r.outfile = c:/bayesx/r-test/pspline/output/reml_gauss_gam4
r.regress norm1 = x1(rw2) + x2(rw1), family=gaussian using d


r.outfile = c:/bayesx/r-test/pspline/output/reml_gauss_vcm1
r.regress norm2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=gaussian using d

r.outfile = c:/bayesx/r-test/pspline/output/reml_gauss_vcm2
r.regress norm2 = x1(psplinerw2) + x2(psplinerw2) + z1 + z1*x3(psplinerw2, center), family=gaussian using d

r.outfile = c:/bayesx/r-test/pspline/output/reml_gauss_vcm3
r.regress norm2 = x1(rw1) + x2(psplinerw2) + z1*x3(rw2), family=gaussian using d


r.outfile = c:/bayesx/r-test/pspline/output/reml_bernoulli_logit_gam1
r.regress bernlogit1 = x1(psplinerw2) + x2(psplinerw2), family=binomial using d

r.outfile = c:/bayesx/r-test/pspline/output/reml_bernoulli_logit_vcm1
r.regress bernlogit2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=binomial using d


r.outfile = c:/bayesx/r-test/pspline/output/reml_bernoulli_probit_gam1
r.regress bernprobit1 = x1(psplinerw2) + x2(psplinerw2), family=binomialprobit using d

r.outfile = c:/bayesx/r-test/pspline/output/reml_bernoulli_probit_vcm1
r.regress bernprobit2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=binomialprobit using d


r.outfile = c:/bayesx/r-test/pspline/output/reml_bernoulli_cloglog_gam1
r.regress berncloglog1 = x1(psplinerw2) + x2(psplinerw2), family=binomialcomploglog using d

r.outfile = c:/bayesx/r-test/pspline/output/reml_bernoulli_cloglog_vcm1
r.regress berncloglog2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=binomialcomploglog using d


r.outfile = c:/bayesx/r-test/pspline/output/reml_binomial_gam1
r.regress bin1 = x1(psplinerw2) + x2(psplinerw2) weight binw, family=binomial using d

r.outfile = c:/bayesx/r-test/pspline/output/reml_binomial_vcm1
r.regress bin2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2) weight binw, family=binomial using d


r.outfile = c:/bayesx/r-test/pspline/output/reml_poisson_gam1
r.regress pois1 = x1(psplinerw2) + x2(psplinerw2), family=poisson using d

r.outfile = c:/bayesx/r-test/pspline/output/reml_poisson_vcm1
r.regress pois2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=poisson using d


% bayesreg

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_gauss_gam1
b.regress norm1 = x1(psplinerw2) + x2(psplinerw2), family=gaussian iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_gauss_gam2
b.regress norm1 = x1(psplinerw2, degree=2, nrknots=12) + x2(psplinerw1, degree=1, nrknots=30), family=gaussian iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_gauss_gam3
b.regress norm1 = x1(psplinerw1, degree=0, nrknots=12) + x2(psplinerw2, degree=1, nrknots=30), family=gaussian iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_gauss_gam4
b.regress norm1 = x1(rw2) + x2(rw1), family=gaussian iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_gauss_gam5
b.regress norm1 = x1(psplinerw2, a=-1, b=0) + x2(psplinerw2, a=0, b=0), family=gaussian iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_gauss_gam6
b.regress norm1 = x1(psplinerw2, a=1, b=0.1) + x2(psplinerw2, a=-0.5, b=0), family=gaussian iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_gauss_gam7
b.regress norm1 = x1(psplinerw2, derivative) + x2(psplinerw2), family=gaussian iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:/bayesx/r-test/pspline/output/mcmc_gauss_vcm1
b.regress norm2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=gaussian iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_gauss_vcm2
b.regress norm2 = x1(psplinerw2) + x2(psplinerw2) + z1 + z1*x3(psplinerw2, center), family=gaussian iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_gauss_vcm3
b.regress norm2 = x1(rw1) + x2(psplinerw2) + z1*x3(rw2), family=gaussian iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:/bayesx/r-test/pspline/output/mcmc_bernoulli_logit_gam1
b.regress bernlogit1 = x1(psplinerw2) + x2(psplinerw2), family=binomial iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_bernoulli_logit_vcm1
b.regress bernlogit2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=binomial iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:/bayesx/r-test/pspline/output/mcmc_bernoulli_probit_gam1
b.regress bernprobit1 = x1(psplinerw2) + x2(psplinerw2), family=binomialprobit iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_bernoulli_probit_vcm1
b.regress bernprobit2 = x1(psplinerw2, lambda=1000) + x2(psplinerw2, lambda=1000) + z1*x3(psplinerw2,lambda=1000), family=binomialprobit iterations=12000 burnin=2000 step=10 setseed=1234 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_binomial_gam1
b.regress bin1 = x1(psplinerw2) + x2(psplinerw2) weight binw, family=binomial iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_binomial_vcm1
b.regress bin2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2) weight binw, family=binomial iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:/bayesx/r-test/pspline/output/mcmc_poisson_gam1
b.regress pois1 = x1(psplinerw2) + x2(psplinerw2), family=poisson iterations=12000 burnin=2000 step=10 setseed=123 using d

b.outfile = c:/bayesx/r-test/pspline/output/mcmc_poisson_vcm1
b.regress pois2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=poisson iterations=12000 burnin=2000 step=10 setseed=123 using d


% stepwisereg

s.outfile = c:/bayesx/r-test/pspline/output/step_gauss_gam1
s.regress norm1 = x1(psplinerw2) + x2(psplinerw2), family=gaussian using d

s.outfile = c:/bayesx/r-test/pspline/output/step_gauss_gam2
s.regress norm1 = x1(psplinerw2, degree=2, nrknots=12) + x2(psplinerw1, degree=1, nrknots=30), family=gaussian using d

s.outfile = c:/bayesx/r-test/pspline/output/step_gauss_gam3
s.regress norm1 = x1(psplinerw1, degree=0, nrknots=12) + x2(psplinerw2, degree=1, nrknots=30), family=gaussian using d

s.outfile = c:/bayesx/r-test/pspline/output/step_gauss_gam4
s.regress norm1 = x1(rw2) + x2(rw1), family=gaussian using d


s.outfile = c:/bayesx/r-test/pspline/output/step_gauss_vcm1
s.regress norm2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=gaussian using d

s.outfile = c:/bayesx/r-test/pspline/output/step_gauss_vcm2
s.regress norm2 = x1(psplinerw2) + x2(psplinerw2) + z1 + z1*x3(psplinerw2, center), family=gaussian using d

s.outfile = c:/bayesx/r-test/pspline/output/step_gauss_vcm3
s.regress norm2 = x1(rw1) + x2(psplinerw2) + z1*x3(rw2), family=gaussian using d


s.outfile = c:/bayesx/r-test/pspline/output/step_bernoulli_logit_gam1
s.regress bernlogit1 = x1(psplinerw2) + x2(psplinerw2), family=binomial using d

s.outfile = c:/bayesx/r-test/pspline/output/step_bernoulli_logit_vcm1
s.regress bernlogit2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=binomial using d


s.outfile = c:/bayesx/r-test/pspline/output/step_bernoulli_probit_gam1
s.regress bernprobit1 = x1(psplinerw2) + x2(psplinerw2), family=binomialprobit using d

s.outfile = c:/bayesx/r-test/pspline/output/step_bernoulli_probit_vcm1
s.regress bernprobit2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=binomialprobit using d


s.outfile = c:/bayesx/r-test/pspline/output/step_binomial_gam1
s.regress bin1 = x1(psplinerw2) + x2(psplinerw2) weight binw, family=binomial using d

s.outfile = c:/bayesx/r-test/pspline/output/step_binomial_vcm1
s.regress bin2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2) weight binw, family=binomial using d


s.outfile = c:/bayesx/r-test/pspline/output/step_poisson_gam1
s.regress pois1 = x1(psplinerw2) + x2(psplinerw2), family=poisson using d

s.outfile = c:/bayesx/r-test/pspline/output/step_poisson_vcm1
s.regress pois2 = x1(psplinerw2) + x2(psplinerw2) + z1*x3(psplinerw2), family=poisson using d