%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading dataset information %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataset d
d.infile, maxobs=5000 using c:\data\zambia.raw
d.describe

d.tabulate sex
d.descriptive bmi

%%%%%%%%%%%%%%%
% Map objects %
%%%%%%%%%%%%%%%

map m
m.infile using c:\data\zambia.bnd
m.describe

m.reorder
m.outfile, replace using c:\data\zambiasort.bnd

m.outfile, replace graph using c:\data\zambiasort.gra
map m1
m1.infile, graph using c:\data\zambiasort.gra
m1.describe
drop m1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bayesian semiparametric regression %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bayesreg b
b.outfile = c:\data\b
logopen, replace using c:\data\logmcmc.txt
b.regress hazstd = rcw + edu1 + edu2 + tpr + sex + bmi(psplinerw2) + agc(psplinerw2) + district(spatial, map=m) + district(random), family=gaussian iterations=12000 burnin=2000 step=10 predict using d
logclose

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualising estimation results %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Post estimation commands

b.plotnonp 1
b.plotnonp 3

b.plotnonp 1, replace outfile = c:\data\f_bmi.ps

b.drawmap 5

% Graph objects

dataset res
res.infile using c:\data\b_f_bmi_pspline.res
graph g
g.plot bmi pmean pqu2p5 pqu10 pqu90 pqu97p5 using res

res.infile using c:\data\b_f_district_spatial.res
g.drawmap pmean district, map=m using res
g.drawmap pcat95 district, map=m using res

res.infile using c:\data\b_f_district_random.res
g.drawmap pmean district, map=m using res

%%%%%%%%%%%%%%%%%%%%%%%%
% Customising graphics %
%%%%%%%%%%%%%%%%%%%%%%%%

b.plotnonp 1, levels=2
b.plotnonp 1, title="Mother body mass index"
b.plotnonp 1, xlab="bmi" ylab="f_bmi" title="Mother body mass index"
b.plotnonp 1, xlab="bmi" ylab="f_bmi" title="Mother body mass index" ylimbottom=-0.8 ylimtop=0.6 ystep=0.2 xlimbottom=12 xlimtop=40

b.plotnonp 3, xlab="age" ylab="f_age" title="Age of the child in months" ylimbottom=-0.3  ystep=0.3 xlimbottom=0 xlimtop=60 xstep=10

b.drawmap 5, color swapcolors

b.drawmap 5, color swapcolors title="Structured spatial effect"
b.drawmap 5, color swapcolors title="Structured spatial effect" lowerlimit=-0.3 upperlimit=0.3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autocorrelation functions and sampling paths %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b.plotautocor, maxlag=250
b.plotautocor, mean

b.getsample

%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity analysis %
%%%%%%%%%%%%%%%%%%%%%%%%

b.regress hazstd = rcw + edu1 + edu2 + tpr + sex + bmi(psplinerw2,a=0.00001,b=0.00001) + agc(psplinerw2,a=0.00001,b=0.00001) + district(spatial, map=m,a=0.00001,b=0.00001) + district(random,a=0.00001,b=0.00001), family=gaussian iterations=12000 burnin=2000 step=10 predict using d
b.plotnonp 1
b.plotnonp 3

b.regress hazstd = rcw + edu1 + edu2 + tpr + sex + bmi(psplinerw2,a=1,b=0.005) + agc(psplinerw2,a=1,b=0.005) + district(spatial,map=m,a=1,b=0.005) + district(random,a=1,b=0.005), family=gaussian iterations=12000 burnin=2000 step=10 predict using d
b.plotnonp 1
b.plotnonp 3

b.regress hazstd = rcw + edu1 + edu2 + tpr + sex + bmi(psplinerw2,a=1,b=0.00005) + agc(psplinerw2,a=1,b=0.00005) + district(spatial,map=m,a=1,b=0.00005) + district(random,a=1,b=0.00005), family=gaussian iterations=12000 burnin=2000 step=10 predict using d
b.plotnonp 1
b.plotnonp 3
