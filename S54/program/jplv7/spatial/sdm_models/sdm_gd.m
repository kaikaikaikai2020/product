% PURPOSE: An example of using sdm_g() Gibbs sampling
%          spatial durbin model
%          (on a small data set)                  
%---------------------------------------------------
% USAGE: sdm_gd (see also sdm_gd2 for a large data set)
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
sige = 2;
k = 2;
x = [ones(n,1) randn(n,k)];
xmat = [x W*x(:,2:end)];
[junk nk] = size(xmat);
beta = ones(nk,1);
y = (IN - rho*W)\(xmat*beta) + (IN - rho*W)\(randn(n,1)*sqrt(sige));

vnames = strvcat('crime','constant','income','hvalue');
info.lflag = 0; % don't use Barry-Pace ln det approximation
res = sdm(y,x,W);
prt(res,vnames);

% Gibbs sampling function homoscedastic prior
prior.rval = 200; % homoscedastic prior for comparison
prior.lflag = 0;
% to maximum likelihood estimates
ndraw = 5000;
nomit = 2500;
prior.dflag = 1; % Metropolis-Hastings sampling for rho
results = sdm_g(y,x,W,ndraw,nomit,prior);
results.tflag = 'tstat';

prt(results,vnames);

prior2.novi = 1;
prior2.lflag = 0;
% uses default numerical integration for rho
results2 = sdm_g(y,x,W,ndraw,nomit,prior2);
results2.tflag = 'tstat';

prt(results2,vnames);


[h1,f1,y1] = pltdens(results.pdraw);
[h2,f2,y2] = pltdens(results2.pdraw);
plot(y1,f1,'.r',y2,f2,'.g');
legend('metropolis-h','integration');
title('posterior distributions for rho');
xlabel('rho values');

