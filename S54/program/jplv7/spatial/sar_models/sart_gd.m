% PURPOSE: An example of using sart_g() Gibbs sampling
%          spatial autoregressive tobit model
%          (on a small data set)                  
%---------------------------------------------------
% USAGE: sart_gd (see also sart_gd2 for a large data set)
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
sige = 0.1;
k = 3;
x = randn(n,k);
beta(1,1) = -1.0;
beta(2,1) = 1.0;
beta(3,1) = 1.0;

y = (IN-rho*W)\(x*beta) + (IN-rho*W)\(randn(n,1)*sqrt(sige)); 
ysave = y;

res = sar(ysave,x,W); % maximum likelihood estimates
fprintf(1,'max like estimates based on actual y-values \n');
prt(res);             % based on non-truncated data

limit = 0;
ind = find(y < limit);
if length(ind) > 0
y(ind,1) = limit; % censor  values
else
    error('no censored values');
end;


%prior.rval = 4;   % heteroscedastic prior
prior.novi = 1;  % homoscedastic prior
ndraw = 2500;
nomit = 500;

result = sar_g(ysave,x,W,ndraw,nomit,prior); % MCMC estimates 
fprintf(1,'Bayesian MCMC estimates based on actual y-values \n');
prt(result);

prior2.novi = 1;  % homoscedastic prior
prior2.limit = limit;
prior2.trunc = 'left';

result2 = sart_g(y,x,W,ndraw,nomit,prior2);
fprintf(1,'Bayesian MCMC estimates based on truncated y-values \n');
prt(result2);

tt=1:n;
plot(tt,ysave,tt,result2.ymean,'--');
title('actual y vs mean of latent y-draws');
pause;

% plot densities for comparison
[h1,f1,y1] = pltdens(result.bdraw(:,1));
[h2,f2,y2] = pltdens(result2.bdraw(:,1));
[h3,f3,y3] = pltdens(result.bdraw(:,2));
[h4,f4,y4] = pltdens(result2.bdraw(:,2));
[h5,f5,y5] = pltdens(result.bdraw(:,3));
[h6,f6,y6] = pltdens(result2.bdraw(:,3));

plot(y1,f1,'.r',y2,f2,'.g');
legend('sar','sart');
xlabel(['true b =' num2str(beta(1,1))]);
pause;
plot(y3,f3,'.r',y4,f4,'.g');
legend('sar','sart');
xlabel(['true b =' num2str(beta(2,1))]);
pause;
plot(y5,f5,'.r',y6,f6,'.g');
legend('sar','sart');
xlabel(['true b =' num2str(beta(3,1))]);
pause;


[h5,f5,y5] = pltdens(result.pdraw);
[h6,f6,y6] = pltdens(result2.pdraw);

plot(y5,f5,'.r',y6,f6,'.g');
legend('sar','sart');
xlabel(['true rho =' num2str(rho)]);



