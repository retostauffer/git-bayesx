clear
set memory 5m

infile region x y using c:\bayesx\test\testdata\kreisecentroid.raw

replace x = (x-2323.031)/1201.466
replace y = (y-3504.746)/1968.201

gen fx = 0.4*(x+y)
sum fx
replace fx = fx-_result(3)

expand 4

gen x2= -1+2*(uniform()<=0.5)


gen eta = 1+x2+fx*x2


gen y1 = eta + 0.3*invnorm(uniform())

outsheet using c:\bayesx\test\testdata\gaussian_spatial_vcm_center.raw  , replace


