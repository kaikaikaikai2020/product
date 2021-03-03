% PURPOSE: An example of using sdmp_g() Gibbs sampling
%          spatial durbin probit model
%          (on a small data set)                  
%---------------------------------------------------
% USAGE: sdmp_gd (see also sdmp_gd2 for a large data set)
%---------------------------------------------------

clear all;

% W-matrix from Anselin's neigbhorhood crime data set
load anselin.dat; % standardized 1st-order spatial weight matrix
latt = anselin(:,4);
long = anselin(:,5);
W = make_neighborsw(latt,long,5); % 5 nearest neighbors weight matrix

[n junk] = size(W);
IN = eye(n); 
rho = 0.8;  % true value of rho
sige = 10;
k = 3;

x = randn(n,k);
beta(1,1) = -0.5;
beta(2,1) = 2.0;
beta(3,1) = -2.0;

y = (IN-rho*W)\(x*beta) + (IN-rho*W)\(randn(n,1)*sqrt(sige)); 


z = (y > 0);
z = ones(n,1).*z; % eliminate a logical vector

% Gibbs sampling function homoscedastic prior
prior.novi = 1; % homoscedastic prior for comparison
% to maximum likelihood estimates based on true y-values
ndraw = 2500;
nomit = 500;

results = sdm_g(y,x,W,ndraw,nomit,prior);
fprintf(1,'results based on actual y-values');
vnames = strvcat('y','x1','x2','x3');
prt(results,vnames);


results2 = sdmp_g(z,x,W,ndraw,nomit,prior);
fprintf(1,'results based on 0,1 y-values');
prt(results2,vnames);


tt=1:n;
plot(tt,y,'-b',tt,results2.yhat,'-r',tt,results.yhat,'-g');
legend('actual','sdmp predicted','sdm predicted');
pause;

[ysort yind] = sort(z);

tt=1:n;
plot(tt,ysort,'or',tt,results2.yprob(yind,1),'.b');
title('probabilities');
legend('actual 0,1 values','probabilities');


% plot densities for comparison
[h1,f1,y1] = pltdens(results.bdraw(:,1));
[h2,f2,y2] = pltdens(results2.bdraw(:,1));
[h3,f3,y3] = pltdens(results.bdraw(:,2));
[h4,f4,y4] = pltdens(results2.bdraw(:,2));
[h5,f5,y5] = pltdens(results.bdraw(:,3));
[h6,f6,y6] = pltdens(results2.bdraw(:,3));

plot(y1,f1,'.r',y2,f2,'.g');
legend('sar','sarp');
xlabel(['true b =' num2str(beta(1,1))]);
pause;
plot(y3,f3,'.r',y4,f4,'.g');
legend('sar','sarp');
xlabel(['true b =' num2str(beta(2,1))]);
pause;
plot(y5,f5,'.r',y6,f6,'.g');
legend('sar','sarp');
xlabel(['true b =' num2str(beta(3,1))]);
pause;


[h5,f5,y5] = pltdens(results.pdraw);
[h6,f6,y6] = pltdens(results2.pdraw);

plot(y5,f5,'.r',y6,f6,'.g');
legend('sar','sarp');
xlabel(['true rho =' num2str(rho)]);
pause;

