
% Reading the data into BayesX

dataset credit
credit.infile using c:\bayes\examples\credit.raw
credit.generate account1  = 1*(account=1)-1*(account=3)
credit.generate account2  = 1*(account=2)-1*(account=3)
credit.generate payment1 = 1*(payment=1)-1*(payment=2)
credit.generate intuse1 = 1*(intuse=1)-1*(intuse=2)
credit.generate marstat1 = 1*(marstat=1)-1*(marstat=2)



% Creating a bayesreg object 

bayesreg b
b.outfile = c:\results\credit



% Probit model: Estimating fixed effects 

b.regress  y = account1 + account2 + duration + amount + payment1 + intuse1 + marstat1  , predict iterations=6000 burnin=1000 step=5 family=binomialprobit using credit



% Probit model: Nonlinear effects for duration and amount 

b.regress  y = account1 + account2 + duration(psplinerw2) + amount(psplinerw2) + payment1 + intuse1 + marstat1  , predict iterations=6000 burnin=1000 step=5 family=binomialprobit using credit

b.plotnonp 1 , outfile="c:\results\credit_duration.ps"
b.plotnonp 3 , outfile="c:\results\credit_amount.ps"

b.plotnonp 1 , outfile="c:\results\credit_duration.ps" replace xlab="duration" ylab="f(duration)" title="effect of duration"
b.plotnonp 3 , outfile="c:\results\credit_amount.ps" replace xlab="amount" ylab="f(amount)" title="effect of amount"

b.plotautocor , outfile="c:\results\credit_autocor.ps"


% Logit model

b.regress  y = account1 + account2 + duration(psplinerw2) + amount(psplinerw2) + payment1 + intuse1 + marstat1  , predict iterations=6000 burnin=1000 step=5 family=binomial using credit
b.plotnonp 1 , outfile="c:\results\credit_logit_duration.ps" replace xlab="duration" ylab="f(duration)" title="effect of duration"
b.plotnonp 3 , outfile="c:\results\credit_logit_amount.ps" replace xlab="amount" ylab="f(amount)" title="effect of amount"

b.plotautocor , outfile="c:\results\credit_logit_autocor.ps"


% Logit model with IWLS proposals based on posterior modes

b.regress  y = account1 + account2 + duration(psplinerw2,proposal=iwlsmode) + amount(psplinerw2,proposal=iwlsmode) + payment1 + intuse1 + marstat1  , predict iterations=6000 burnin=1000 step=5 family=binomial using credit


% Probit model: Varying the hyperparameters

% a=0.0001 b=0.0001
b.regress  y = account1 + account2 + duration(psplinerw2,a=0.0001,b=0.0001) + amount(psplinerw2,a=0.0001,b=0.0001) + payment1 + intuse1 + marstat1  , predict iterations=6000 burnin=1000 step=5 family=binomialprobit using credit
b.plotnonp 1 , outfile="c:\results\credit_duration_a0001b0001.ps" replace xlab="duration" ylab="f(duration)" title="a: 0.0001 and b: 0.0001"
b.plotnonp 3 , outfile="c:\results\credit_amount_a0001b0001.ps" replace xlab="amount" ylab="f(amount)" title="a: 0.0001 and b: 0.0001"

% a=1.0 b=0.005
b.regress  y = account1 + account2 + duration(psplinerw2,a=1.0,b=0.005) + amount(psplinerw2,a=1.0,b=0.005) + payment1 + intuse1 + marstat1  , predict iterations=6000 burnin=1000 step=5 family=binomialprobit using credit
b.plotnonp 1 , outfile="c:\results\credit_duration_a1b005.ps" replace xlab="duration" ylab="f(duration)" title="a: 1.0 and b: 0.005"
b.plotnonp 3 , outfile="c:\results\credit_amount_a1b005.ps" replace xlab="amount" ylab="f(amount)" title="a: 1.0 and b: 0.005"


