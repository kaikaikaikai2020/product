% PURPOSE: An example of using sarp_g() Gibbs sampling
%          spatial autoregressive probit model
%          (on a small data set)                  
%---------------------------------------------------
% USAGE: sarp_gd (see also sarp_gd2 for a large data set)
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

disp('maximum likelihood estimates based on continuous y');
res = sar(y,x,W);
prt(res);

z = (y > 0);
z = ones(n,1).*z; % eliminate a logical vector

% Gibbs sampling function homoscedastic prior
% to maximum likelihood estimates
ndraw = 6000;
nomit = 1000;
prior.novi = 1; % homoscedastic prior
prior.rmin = -1;
prior.rmax = 1;

result = sar_g(y,x,W,ndraw,nomit,prior);

%prior.rval = 4;% heteroscedastic prior
result2 = sarp_g(z,x,W,ndraw,nomit,prior);

% plot draws for comparison
tt=1:ndraw-nomit;
plot(tt,result.bdraw(:,1),'-r',tt,result2.bdraw(:,1),'-g');
xlabel(['true b =' num2str(beta(1,1))]);
legend('sar','sarp');
pause;
plot(tt,result.bdraw(:,2),'-r',tt,result2.bdraw(:,2),'-g');
xlabel(['true b =' num2str(beta(2,1))]);
legend('sar','sarp');
pause;
plot(tt,result.bdraw(:,3),'-r',tt,result2.bdraw(:,3),'-g');
xlabel(['true b =' num2str(beta(3,1))]);
legend('sar','sarp');
pause;

plot(tt,result.pdraw,'-r',tt,result2.pdraw,'-g');
xlabel(['true rho =' num2str(rho)]);
legend('sar','sarp');
pause;


% plot densities for comparison
[h1,f1,y1] = pltdens(result.bdraw(:,1));
[h2,f2,y2] = pltdens(result2.bdraw(:,1));
[h3,f3,y3] = pltdens(result.bdraw(:,2));
[h4,f4,y4] = pltdens(result2.bdraw(:,2));
[h5,f5,y5] = pltdens(result.bdraw(:,3));
[h6,f6,y6] = pltdens(result2.bdraw(:,3));

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


[h5,f5,y5] = pltdens(result.pdraw);
[h6,f6,y6] = pltdens(result2.pdraw);

plot(y5,f5,'.r',y6,f6,'.g');
legend('sar','sarp');
xlabel(['true rho =' num2str(rho)]);
pause;

tt=1:n;
plot(tt,y,'-r',tt,result.yhat,'-g',tt,result2.ymean,'-b');
legend('actual y','sar yhat','sarp mean of y draws');
pause;

plot(tt,z,'o',tt,result2.yprob,'.');
legend('0-1 y-values,','probabilities');


prt(result);
prt(result2);
