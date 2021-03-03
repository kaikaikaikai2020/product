% PURPOSE: An example of using sac() on a small dataset
%          general spatial model
%                              
%---------------------------------------------------
% USAGE: sac_d
%---------------------------------------------------

clear all;

% load Anselin (1988) Columbus neighborhood crime data
load anselin.dat; 
n = length(anselin);
randn('seed',202020);
x = randn(n,3);
xc = anselin(:,4);
yc = anselin(:,5);
vnames = strvcat('crime','constant','income','hvalue');
[j1 W j2] = xy2cont(xc,yc);

% create a nearest neighbor weight matrix for the error spatial autocorrelation
W2 = make_nnw(xc,yc,2);

% do Monte Carlo generation of an SAC model
sige = 10; 
evec = randn(n,1)*sqrt(sige);
beta = ones(3,1);
rho = 0.4; lam = -0.7; 
A = eye(n) - rho*W;  
B = eye(n) - lam*W2; 
u = B\evec;
y = A\(x*beta) + A\u; % generate some data

res = sac(y,x,W,W2);
% print the output with variable names
prt(res,vnames);
plt(res);

info.lflag = 0;
res2 = sar(y,x,W,info);
prt(res2,vnames);

res3 = sem(y,x,W,info);
prt(res3,vnames);


