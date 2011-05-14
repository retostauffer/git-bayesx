% usefile c:\bayesx\r-test\linear\r-test.prg

dataset d
d.infile using c:\bayesx\r-test\linear\data.raw

remlreg r
bayesreg b
stepwisereg s


% remlreg

r.outfile = c:\bayesx\r-test\linear\output\reml_gauss_linear1
r.regress norm1 = x1 + x2, family=gaussian using d


r.outfile = c:\bayesx\r-test\linear\output\reml_bernoulli_logit_linear1
r.regress bernlogit1 = x1 + x2, family=binomial using d


r.outfile = c:\bayesx\r-test\linear\output\reml_bernoulli_probit_linear1
r.regress bernprobit1 = x1 + x2, family=binomialprobit using d


r.outfile = c:\bayesx\r-test\linear\output\reml_bernoulli_cloglog_linear1
r.regress berncloglog1 = x1 + x2, family=binomialcomploglog using d


r.outfile = c:\bayesx\r-test\linear\output\reml_binomial_linear1
r.regress bin1 = x1 + x2 weight binw, family=binomial using d


r.outfile = c:\bayesx\r-test\linear\output\reml_poisson_linear1
r.regress pois1 = x1 + x2, family=poisson using d


r.outfile = c:\bayesx\r-test\linear\output\reml_cumprobit_linear1
r.regress cump = x1 + x2, family=cumprobit using d


r.outfile = c:\bayesx\r-test\linear\output\reml_multlogit_linear1
r.regress mlog = x1 + x2, family=multinomial reference=3 using d



% bayesreg

b.outfile = c:\bayesx\r-test\linear\output\mcmc_gauss_linear1
b.regress norm1 = x1 + x2, family=gaussian iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:\bayesx\r-test\linear\output\mcmc_bernoulli_logit_linear1
b.regress bernlogit1 = x1 + x2, family=binomial iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:\bayesx\r-test\linear\output\mcmc_bernoulli_probit_linear1
b.regress bernprobit1 = x1 + x2, family=binomialprobit iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:\bayesx\r-test\linear\output\mcmc_binomial_linear1
b.regress bin1 = x1 + x2 weight binw, family=binomial iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:\bayesx\r-test\linear\output\mcmc_poisson_linear1
b.regress pois1 = x1 + x2, family=poisson iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:\bayesx\r-test\linear\output\mcmc_cumprobit_linear1
b.regress cump = x1 + x2, family=cumprobit iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:\bayesx\r-test\linear\output\mcmc_multlogit_linear1
b.regress mlog = x1 + x2, family=multinomial reference=3 iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:\bayesx\r-test\linear\output\mcmc_multprobit_linear1
b.regress mlog = x1 + x2, family=multinomialprobit reference=3 iterations=12000 burnin=2000 step=10 setseed=123 using d


b.outfile = c:\bayesx\r-test\linear\output\mcmc_negbinomial_linear1
b.regress  nbin = x1 + x2, family=nbinomial distopt=nb iterations=12000 burnin=2000 step=10 setseed=123 using d


% b.outfile = c:\bayesx\r-test\linear\output\mcmc_zip_linear1
% b.regress zip = x1 + x2, family=zip zipdistopt=zip iterations=12000 burnin=2000 step=10 setseed=123 using d



% stepwisereg

s.outfile = c:\bayesx\r-test\linear\output\step_gauss_linear1
s.regress norm1 = x1 + x2, family=gaussian CI=MCMCbootstrap step=10 setseed=123 using d


s.outfile = c:\bayesx\r-test\linear\output\step_bernoulli_logit_linear1
s.regress bernlogit1 = x1 + x2, family=binomial CI=MCMCbootstrap step=10 setseed=123 using d


% s.outfile = c:\bayesx\r-test\linear\output\step_bernoulli_probit_linear1
% s.regress bernprobit1 = x1 + x2, family=binomialprobit CI=MCMCbootstrap step=10 setseed=123 using d


s.outfile = c:\bayesx\r-test\linear\output\step_binomial_linear1
s.regress bin1 = x1 + x2 weight binw, family=binomial CI=MCMCbootstrap step=10 setseed=123 using d


s.outfile = c:\bayesx\r-test\linear\output\step_poisson_linear1
s.regress pois1 = x1 + x2, family=poisson CI=MCMCbootstrap step=10 setseed=123 using d


s.outfile = c:\bayesx\r-test\linear\output\step_multlogit_linear1
s.regress mlog = x1 + x2, family=multinomial reference=3 using d



