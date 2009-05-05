
clear
set obs 1000

gen x1=-3+6*uniform()
gen x2=-3+6*uniform()

gen eta = -0.25+0.4*x1-0.6*x2
gen pi = exp(eta)/(1+exp(eta))
gen y = 1*(uniform() <= pi)


outsheet using c:\bayesx\testh\testdata\hbinomial_linear.raw , replace
save c:\bayesx\testh\testdata\hbinomial_linear.dta , replace






