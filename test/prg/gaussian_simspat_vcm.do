clear
set memory 5m

infile region x y using c:\bayesx\test\testdata\kreisecentroid.raw

replace x = (x-2323.031)/1201.466
replace y = (y-3504.746)/1968.201

gen x2 = uniform()

gen fx = 0.4*(x+y)*x2
sum fx
replace fx = fx-_result(3)

gen eta = 1+x2+fx

expand 4

gen y1 = eta + 0.3*invnorm(uniform())

outsheet using c:\bayesx\test\testdata\gaussian_spatial_vcm_center.raw  , replace


