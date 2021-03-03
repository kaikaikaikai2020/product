% PURPOSE: An example of using sem_g()
%          Gibbs sampling spatial autoregressive model
%          on a large data set                    
%---------------------------------------------------
% USAGE: sem_gd2 (see sem_gd for a small data set)
%---------------------------------------------------

clear all;
% NOTE a large data set with 3107 observations
% from Pace and Barry, takes around 150-250 seconds
load elect.dat;             % load data on votes
y =  log(elect(:,7)./elect(:,8));
x1 = log(elect(:,9)./elect(:,8));
x2 = log(elect(:,10)./elect(:,8));
x3 = log(elect(:,11)./elect(:,8));
n = length(y); 
x = [ones(n,1) x1 x2 x3];
clear x1; clear x2; clear x3;
xc = elect(:,5);
yc = elect(:,6);
[j1 W j2] = xy2cont(xc,yc);
clear elect;                % conserve on RAM memory
n = 3107;
vnames = strvcat('voters','const','educ','homeowners','income');
res = sem(y,x,W);
prt(res,vnames);

% do Gibbs sampling estimation
ndraw = 2500; 
nomit = 500;
prior.novi = 1; % homoscedastic prior
% uses default M-H sampling for rho
% uses default Pace-Barry lndet approximation

resg = sem_g(y,x,W,ndraw,nomit,prior);
resg.tflag = 'tstat';
% these homoscedastic results should match max lik results
prt(resg,vnames);

prior2.rval = 4; % heteroscedastic prior
% uses default M-H sampling for rho
% uses default Pace-Barry lndet approximation

resg2 = sem_g(y,x,W,ndraw,nomit,prior2);
resg2.tflag = 'tstat'; % print bogus t-statistics for comparison
prt(resg2,vnames);     % with maximum likelihood estimates


[h1,f1,y1] = pltdens(resg.pdraw);
[h2,f2,y2] = pltdens(resg2.pdraw);
plot(y1,f1,'.r',y2,f2,'.g');
legend('homoscedastic','heteroscedastic');
title('rho posterior distributions');

probs = model_probs(resg,resg2);

probs

