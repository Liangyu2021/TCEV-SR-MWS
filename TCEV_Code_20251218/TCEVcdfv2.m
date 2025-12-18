function ycalc=TCEVcdfv2(x,par)
% It returns the CDF according to TCEV function in x according to 
% the 4 parameters included in 'par'
% x can be a vector
% par=(alfa1, beta1 (1ªdistribu) and alfa2,beta2 (2ª %distribu))

% par must be a vector (1x4)

alfa1=par(1);
beta1=par(2);
alfa2=par(3);
beta2=par(4);

ycalc=exp(-exp(-alfa1*(x-beta1))).*exp(-exp(-alfa2*(x-beta2)));
