## Namibian hake data from Polacheck et al., 1993
## Fitting Surplus Production Models: Comparing Methods and Measuring Uncertainty
years=1965:1988

### catch data (MT)
Ct=c(93.510, 212.444, 195.032, 382.712, 320.430, 402.467, 365.557, 606.084, 377.642,
     318.836, 309.374, 389.020, 276.901, 254.251, 170.006, 97.181, 90.523, 176.532, 216.181, 228.672, 212.177,
     231.179, 136.942, 212.000)

## cpue (MT/standardised trawl hour)
It=c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55,
     0.53, 0.58, 0.64, 0.66, 0.65, 0.61, 0.63)


ntimes=length(years)

names(Ct)<-years
names(It)<-years



######## biological info (for ASPM) from Forrest et al., 2007  ###########

########################
###### Input pars ######
########################

a_mat=4
a_sel=3
nages=12

sig_sel=0.2
#sig_sel=0.1*a_sel
sig_mat=0.2*a_mat

Linf=111
kappa=0.14
a0=0
lwa=0.00001 *1e-6
lwb=3
#sigR=0.4
sigR=0.3852532
FemaleProp=0.5

## assumed median values of M and steepness parameters
#M=1.5*kappa
M=0.35
h=0.8
