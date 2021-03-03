% PURPOSE: An example of using sdmp_g() on a large data set   
%          Gibbs sampling spatial durbin probit model                         
%---------------------------------------------------
% USAGE: sdmp_gd2 (see sdmp_gd for a small data set)
%---------------------------------------------------

clear all;
% NOTE a large data set with 3107 observations
% from Pace and Barry, 
load elect.dat;             % load data on votes
latt = elect(:,5);
long = elect(:,6);
n = length(latt);
k = 2;
x = [randn(n,k)];
clear elect;                % conserve on RAM memory
n = 3107;
[junk W junk] = xy2cont(latt,long);

rho = 0.7;
beta(1,1) = 1.0;
beta(2,1) = -1.0;

sige = 1;

y = (speye(n) - rho*W)\(x*beta) + (speye(n) - rho*W)\randn(n,1)*sqrt(sige);
ysave = y;
y = (y > mean(y));

% do Gibbs sampling estimation
ndraw = 1200; 
nomit = 200;
%prior.rval = 200;
prior.novi = 1;

res = sdm_gc(ysave,x,W,ndraw,nomit,prior); % MCMC estimates based on
prt(res);                                  % non-truncated data for comparison

resg = sdmp_g(y,x,W,ndraw,nomit,prior); % tobit estimates
prt(resg);

tt=1:n;
plot(tt,ysave,'.b',tt,resg.yhat,'.r',tt,res.yhat,'.g');
legend('actual','sdmp predicted','sdm predicted');
pause;

[ysort yind] = sort(y);

tt=1:n;
plot(tt,ysort,'or',tt,resg.yprob(yind,1),'.g');