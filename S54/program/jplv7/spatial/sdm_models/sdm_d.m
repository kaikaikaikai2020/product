% PURPOSE: An example of using sdm() max likelihood
%          estimation of the spatial durbin model
%          (on a small data set)                  
%---------------------------------------------------
% USAGE: sdm_d (see also sdm_d2 for a large data set)
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
k = 3;
x = [ones(n,1) anselin(:,2:3)];
y = anselin(:,1);
vnames = strvcat('crime','constant','income','hvalue');
info.lflag = 0;
info.rmin = -1;
info.rmax = 1;
res = sdm(y,x,W,info);
prt(res,vnames);

% generate our own model without an intercept term
x = [randn(n,k)];
xsdm = [x W*x];
[n,nvar] = size(xsdm);
beta = 2*ones(nvar,1);
for i=1:2:nvar
beta(i,1) = -2.0;
end;

vnames = strvcat('y','x1','x2','x3');

y = (IN-rho*W)\(xsdm*beta) + (IN-rho*W)\(randn(n,1)*sqrt(sige)); 

res2 = sdm(y,x,W,info);
prt(res2,vnames);
