
clear
set obs 1000

generate x1 = -3+6*uniform()
generate x2 = -3+6*uniform()

generate eta1 = 0.5*x1-0.5*x2
generate eta2 = 0.5*x1+0.5*x2
gen u0 = 0
gen u1 = eta1+invnorm(uniform())
gen u2 = eta2+invnorm(uniform())
gen y1 = 1*(u1 > u2 & u1 > u0) 
gen y2 = 1*(u2 > u1 & u2 > u0) 

gen y = 1*(u1 > u2 & u1 > u0)+2*(u2 > u1 & u2 > u0)

mprobit y x1 x2 , baseoutcome(0)

outsheet using c:\bayesx\testh\testdata\multinomprobit.raw

 