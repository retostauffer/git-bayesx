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
% Stepwise semiparametric regression %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delimiter = ;

stepwisereg s;


s.outfile = c:\data\s;
logopen, replace using c:\data\logstep.txt;

s.regress hazstd = rcw + edu(factor,reference=0,coding=effect) + tpr + sex + 
bmi(psplinerw2) + agc(psplinerw2) + 
district(spatial,map=m) + district(random), family=gaussian predict using d;

logclose;


s.outfile = c:\bayesx\output\sc;

s.regress hazstd = rcw + edu(factor,reference=0,coding=effect) + tpr + sex + 
bmi(psplinerw2) + agc(psplinerw2) + 
district(spatial, map=m) + district(random), 
CI=MCMCselect step=10 iterations=10000 family=gaussian predict using d;


s.outfile = c:\bayesx\output\su;

s.regress hazstd = rcw + edu(factor,reference=0,coding=effect) + tpr + sex + 
bmi(psplinerw2) + agc(psplinerw2) + 
district(spatial, map=m) + district(random), 
CI=MCMCbootstrap bootstrapsamples=99 step=10 iterations=10000 family=gaussian predict using d;

delimiter=return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualising estimation results %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Post estimation commands

s.plotnonp 1

s.plotnonp 1, replace outfile = c:\data\f_bmi.ps


s.plotnonp 2

s.plotnonp 2, replace outfile = c:\data\f_age.ps


s.drawmap 3

% Graph objects

dataset res
res.infile using c:\data\b_f_age_pspline.res
graph g
g.plot age pmean pqu2p5 pqu10 pqu90 pqu97p5 using res

res.infile using c:\data\b_f_district_spatial.res
g.drawmap pmean district, map=m using res
g.drawmap pcat95 district, map=m using res

res.infile using c:\data\b_f_district_random.res
g.drawmap pmean district, map=m using res

%%%%%%%%%%%%%%%%%%%%%%%%
% Customising graphics %
%%%%%%%%%%%%%%%%%%%%%%%%

s.plotnonp 1, title="Mother's body mass index"
s.plotnonp 1, xlab="bmi" ylab="f_bmi" title="Mother's body mass index"
s.plotnonp 1, xlab="bmi" ylab="f_bmi" title="Mother's body mass index" ylimbottom=-0.8 ylimtop=0.6 ystep=0.2 xlimbottom=12 xlimtop=40

s.plotnonp 2, levels=2
s.plotnonp 2, title="Age of the child in months"
s.plotnonp 2, xlab="age" ylab="f_age" title="Age of the child in months"
s.plotnonp 2, xlab="age" ylab="f_age" title="Age of the child in months" ylimbottom=-0.4  ystep=0.4 xlimbottom=0 xlimtop=60 xstep=10

s.drawmap 3, color swapcolors
s.drawmap 3, color swapcolors title="Structured spatial effect"
s.drawmap 3, color swapcolors title="Structured spatial effect" lowerlimit=-0.3 upperlimit=0.3

