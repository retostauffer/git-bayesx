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

m.outfile, replace graph using c:\data\zambia.gra
map m1
m1.infile, graph using c:\data\zambia.gra
m1.describe
drop m1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bayesian semiparametric regression %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

remlreg r
r.outfile = c:\data\r
logopen, replace using c:\data\logreml.txt
r.regress hazstd = rcw + edu1 + edu2 + tpr + sex + bmi(psplinerw2) + agc(psplinerw2) + district(spatial, map=m) + district(random), family=gaussian lowerlim=0.01 eps=0.0005 using d
logclose

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualising estimation results %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Post estimation commands

r.plotnonp 1
r.plotnonp 2

r.plotnonp 1, replace outfile = c:\data\f_bmi.ps

r.drawmap 3

% Graph objects

dataset res
res.infile using c:\data\r_f_bmi_pspline.res
graph g
g.plot bmi pmode ci95lower ci80lower ci80upper ci95upper using res

res.infile using c:\data\r_f_district_spatial.res
g.drawmap pmode district, map=m using res
g.drawmap pcat95 district, map=m using res

res.infile using c:\data\r_f_district_random.res
g.drawmap pmode district, map=m using res

%%%%%%%%%%%%%%%%%%%%%%%%
% Customising graphics %
%%%%%%%%%%%%%%%%%%%%%%%%

r.plotnonp 1, levels=2
r.plotnonp 1, title="Mother body mass index"
r.plotnonp 1, xlab="bmi" ylab="f_bmi" title="Mother body mass index"
r.plotnonp 1, xlab="bmi" ylab="f_bmi" title="Mother body mass index" ylimbottom=-0.3 ylimtop=0.6 ystep=0.3 xlimbottom=12 xlimtop=40

r.plotnonp 2, xlab="age" ylab="f_age" title="Age of the child in months" ylimbottom=-0.3  ystep=0.3 xlimbottom=0 xlimtop=60 xstep=10

r.drawmap 3, color swapcolors
r.drawmap 3, color swapcolors title="Structured spatial effect"
r.drawmap 3, color swapcolors title="Structured spatial effect" lowerlimit=-0.3 upperlimit=0.3

