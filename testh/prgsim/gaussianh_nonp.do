
clear
set obs 100

gen x1=-3+6*(_n-1)/99

gen u = uniform()
sort u
gen x2=-3+6*(_n-1)/99

gen id = _n

gen f1 = x1
gen f2 = sin(x2)
sum f2 
replace f2 = f2-_result(3)
gen etare = f1+f2

gen re = 0.66*invnorm(uniform())
sum re
replace re = re -_result(3)

outsheet using c:\bayesx\testh\testdata\hgaussian_nonp1_re.raw , replace
save c:\bayesx\testh\testdata\hgaussian_nonp1_re.dta , replace

expand 10

gen x3 = round(-3+6*uniform(),0.005) 
gen x4 = round(6*uniform(),0.005) 

gen f3 = cos(x3)
sum f3
replace f3 = f3- _result(3)

gen f4 = log(x4)
sum f4
replace f4 =  f4 - _result(3)

gen eta = 1+f1+f2+f3+f4+re


local i=1
while `i' < 251 {
generate y`i' = eta+0.75*invnorm(uniform())
local i=`i'+1 
}


outsheet using c:\bayesx\testh\testdata\hgaussian_nonp1.raw , replace
save c:\bayesx\testh\testdata\hgaussian_nonp1.dta , replace






