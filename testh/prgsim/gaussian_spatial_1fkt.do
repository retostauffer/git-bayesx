clear
set memory 5m

infile region x y using c:\bayesx\testh\testdata\kreisecentroid.raw

replace x = (x-2323.031)/1201.466
replace y = (y-3504.746)/1968.201

gen fx = 0.4*(x+y)

sum fx

outsheet using c:\bayesx\testh\testdata\gaussian_spatial_1fkt_truefunction.raw  , replace
save c:\bayesx\testh\testdata\gaussian_spatial_1fkt_truefunction.dta  , replace

expand 4

gen y1 = fx + 0.3*invnorm(uniform())

outsheet using c:\bayesx\testh\testdata\gaussian_spatial_1fkt.raw  , replace
save c:\bayesx\testh\testdata\gaussian_spatial_1fkt.dta  , replace



