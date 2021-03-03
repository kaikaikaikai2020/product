% PURPOSE: An example of using sac_gmm
% GM estimation of general spatial model (on a small data set)  
%                                   
%---------------------------------------------------
% USAGE: sac_gmmd (see also sac_gmmd2 for a large data set)
%---------------------------------------------------

clear all;
% W-matrix from Anselin's neigbhorhood crime data set
load anselin.dat; 
xc = anselin(:,4);
yc = anselin(:,5);
% crate standardized 1st-order spatial weight matrix
W = make_neighborsw(xc,yc,3);
[j,M,j] = xy2cont(xc,yc);
[n junk] = size(W);
IN = eye(n); 
rho = 0.6;  % true value of rho
lam = -0.4;  % true value of lam
sige = 1;
k = 4;
x = randn(n,k);
beta = zeros(k,1);
beta(1,1) = 1.0;
beta(2,1) = -1.0;
beta(3,1) = 1.0;
beta(4,1) = -1.0;

u = (IN - lam*M)\(randn(n,1)*sqrt(sige));
B = (IN - rho*W);
y = B\(x*beta) + B\u;

results0 = sac(y,x,W,W);
prt(results0);

results1 = sac_gmm(y,x,W);
prt_gmm(results1);
