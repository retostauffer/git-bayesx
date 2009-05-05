
clear
set obs 1000

gen x1=6*uniform()
gen x2=-3+6*uniform()

gen eta = 2+0.8*x1-0.9*x2

gen y = eta + 0.7*invnorm(uniform())

outsheet using c:\bayesx\testh\testdata\hgaussian_linear.raw , replace
save c:\bayesx\testh\testdata\hgaussian_linear.dta , replace






