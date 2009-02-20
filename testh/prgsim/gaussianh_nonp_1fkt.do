
clear
set obs 1000

gen x1=-3+6*uniform()
gen f1 = sin(x1)

gen y = f1 + 0.7*invnorm(uniform())

outsheet using c:\bayesx\testh\testdata\hgaussian_nonp_1fkt.raw , replace
save c:\bayesx\testh\testdata\hgaussian_nonp_1fkt.dta , replace






