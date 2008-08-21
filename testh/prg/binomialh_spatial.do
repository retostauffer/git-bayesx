clear
set memory 5m

infile region x y using c:\bayesx\testh\testdata\kreisecentroid.raw

replace x = (x-2323.031)/1201.466
replace y = (y-3504.746)/1968.201

gen re = 0.6*invnorm(uniform())
sum re
replace re = re-_result(3)

gen fx = 0.4*(x+y)

sum fx

outsheet using c:\bayesx\testh\testdata\binomialh_spatial_cl.raw  , replace
save c:\bayesx\testh\testdata\binomialh_spatial_cl.dta , replace


expand 4
gen x1=-3+6*uniform()

gen y1 = fx + re + 0.5*invnorm(uniform())


outsheet using c:\bayesx\testh\testdata\gaussianh_spatial.raw  , replace
save c:\bayesx\testh\testdata\gaussianh_spatial.dta , replace



