% PURPOSE: An example of using sar_g() Gibbs sampling
%          spatial autoregressive model
%          (on a small data set)                  
%---------------------------------------------------
% USAGE: sar_gd (see also sar_gd2 for a large data set)
%---------------------------------------------------

clear all;

% W-matrix from Anselin's neigbhorhood crime data set
load anselin.dat; % standardized 1st-order spatial weight matrix
latt = anselin(:,4);
long = anselin(:,5);
[junk W junk] = xy2cont(latt,long);
[n junk] = size(W);
IN = eye(n); 
rho = 0.7;  % true value of rho
sige = 0.5;
k = 3;

x = randn(n,k);
beta(1,1) = -1.0;
beta(2,1) = 1.0;
beta(3,1) = 1.0;

y = (IN-rho*W)\(x*beta) + (IN-rho*W)\(randn(n,1)*sqrt(sige)); 

info.lflag = 0; % don't use Pace-Barry lndet approximation
result0 = sar(y,x,W,info);
prt(result0);

ndraw = 2500;
nomit = 500;
prior.novi = 1; % homoscedastic prior
prior.lflag = 0; % don't use Pace-Barry lndet approximation
results = sar_g(y,x,W,ndraw,nomit,prior);
results.tflag = 'tstat';
prt(results);

prior2.lflag = 0; % don't use Pace-Barry lndet approximation
prior2.rval = 4;  % heteroscedastic model 
results2 = sar_g(y,x,W,ndraw,nomit,prior2);
results2.tflag = 'tstat';
prt(results2);

b1mean = mean(results.bdraw)';
b2mean = mean(results2.bdraw)';

out = [beta result0.beta b1mean b2mean
       sige result0.sige results.sige results2.sige
       rho result0.rho  results.rho results2.rho
       0 result0.rsqr results.rsqr results2.rsqr];
   
in.cnames = strvcat('True Values','Max Lik','Bayes homo','Bayes hetero');
in.rnames = strvcat('Parameters','beta0','beta1','beta2','sigma2','rho','r-squared');

fprintf(1,'\n comparison of estimates \n');
mprint(out,in);

% do a plot of posteriors for rho
[h1,f1,y1] = pltdens(results.pdraw);
[h2,f2,y2] = pltdens(results2.pdraw);

plot(y1,f1,'.r',y2,f2,'.g');
legend('rho posterior homo','rho posterior hetero');
title('rho posterior distributions');


