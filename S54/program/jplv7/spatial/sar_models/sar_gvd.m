% PURPOSE: An example of using sar_gv() Gibbs sampling
%          spatial autoregressive model that constructs
%          a posterior distribution for the heteroscedasticity
%          parameter r (rather than use a degenerate prior on r)
%          (on a small data set)                  
%---------------------------------------------------
% USAGE: sar_gvd 
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


ndraw = 2500;
nomit = 500;
prior.delta = 20; % homoscedastic prior
results = sar_gv(y,x,W,ndraw,nomit,prior);

% construct posterior distribution for r-value from
% sampled values to draw inference about heteroscedasticity

subplot(2,1,1),
hist(results.rdraw);
subplot(2,1,2),
plot(results.acc);
title('acceptance rate for M-H sampling');
fprintf(1,'prior mean for r     = %16.8f \n',results.delta);
fprintf(1,'posterior mean for r = %16.8f \n',mean(results.rdraw));
fprintf(1,'posterior std for r  = %16.8f \n',std(results.rdraw));
pause;

prior2.delta = 4;  % heteroscedastic prior 
results2 = sar_gv(y,x,W,ndraw,nomit,prior2);
subplot(2,1,1),
hist(results2.rdraw);
subplot(2,1,2),
plot(results2.acc);
title('acceptance rate for M-H sampling');
fprintf(1,'prior mean for r     = %16.8f \n',results2.delta);
fprintf(1,'posterior mean for r = %16.8f \n',mean(results2.rdraw));
fprintf(1,'posterior std for r  = %16.8f \n',std(results2.rdraw));
pause;


% do a plot of posteriors for r-value
[h1,f1,y1] = pltdens(results.rdraw);
[h2,f2,y2] = pltdens(results2.rdraw);

subplot(1,1,1);
plot(y1,f1,'.r',y2,f2,'.g');
legend('r-value posterior homo prior','r-value posterior hetero prior');
title('r-parameter posterior distributions');


