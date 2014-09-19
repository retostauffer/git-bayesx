
clear
set obs 1000

gen x1=-3+6*uniform()

gen f1 = 0.5*x
gen y = f1 + 0.7*invnorm(uniform())

outsheet using c:\bayesx\trunk\testh\testdata\hgaussian_nonp_1fkt.raw , replace
save c:\bayesx\trunk\testh\testdata\hgaussian_nonp_1fkt.dta , replace



gen f1 = sin(x1)







