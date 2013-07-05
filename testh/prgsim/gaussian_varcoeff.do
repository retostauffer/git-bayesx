clear
set obs 1000
generate x1 = -3+6*uniform()
generate x2 = 1*(_n > 400)
generate f = sin(x1)
generate fg = x2*f
generate y = fg + 0.3*invnorm(uniform())
outsheet using c:\bayesx\trunk\testh\testdata\gaussian_varcoeff.raw , replace
save c:\bayesx\trunk\testh\testdata\gaussian_varcoeff.dta, replace
