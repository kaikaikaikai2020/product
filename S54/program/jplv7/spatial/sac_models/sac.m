function results = sac(y,x,W1,W2,info)
% PURPOSE: computes general Spatial Model estimates
%  model: y = rho*W1*y + X*b + u,  u = lam*W2*u + e
% ---------------------------------------------------
%  USAGE: results = sac(y,x,W1,W2,info)
%  where: y  = dependent variable vector
%         x  = independent variables matrix
%         W1 = spatial weight matrix (standardized)
%         W2 = spatial weight matrix 
%      info        = an (optional) structure variable with input options
%      info.parm   = (optional) 2x1 vector of starting values for rho, lambda
%      info.convg  = (optional) convergence criterion (default = 1e-4)
%      info.maxit  = (optional) maximum # of iterations (default = 500)
%      info.lmin   = (optional) minimum lambda to search (default = -0.99)
%      info.lmax   = (optional) maximum lambda to search (default = 0.99)
%      info.rmin   = (optional) minimum rho to search (default = -0.99)
%      info.rmax   = (optional) maximum rho to search (default = 0.99)
%      info.lflag  = 0 for full computation (default = 1, fastest)
%                  = 1 for Pace and Barry 1999 MC approximation (fast for very large problems)
%                  = 2 for Pace and Barry 1998 Spline approximation (medium speed)
%      info.order  = order to use with info.lflag = 1 option (default = 50)
%      info.iter   = iterations to use with info.lflag = 1 option (default = 30)     
%      info.hessian = 1 for numerical hessian calculation of var-cov matrix
%                     default = 0 if n < 500, analytical hessian
%                             = 1 if n > 500, numerical hessian
% ---------------------------------------------------
%  RETURNS: a structure 
%         results.meth  = 'sac'
%         results.beta  = bhat
%         results.rho   = rho
%         results.lam   = lambda
%         results.tstat = asymptotic t-stats (last 2 are rho,lambda)
%         results.yhat  = yhat
%         results.resid = residuals
%         results.sige  = sige = e'(I-L*W)'*(I-L*W)*e/n
%         results.rsqr  = rsquared
%         results.rbar  = rbar-squared
%         results.lik   = likelihood function value
%         results.nobs  = nobs
%         results.nvar  = nvars
%         results.y     = y data vector
%         results.iter  = # of iterations taken
%         results.lflag = lflag from input
%         results.liter = info.iter option from input
%         results.order = info.order option from input
%         results.limit = matrix of [rho lower95,logdet approx, upper95] intervals
%                         for the case of lflag = 1
%         results.time1 = time for log determinant calcluation
%         results.time2 = time for eigenvalue calculation
%         results.time3 = time for hessian or information matrix calculation
%         results.time4 = time for optimization
%  --------------------------------------------------
%  SEE ALSO: prt_spat(results), prt
% ---------------------------------------------------
% REFERENCES: Luc Anselin Spatial Econometrics (1988) 
%             pages 64-65 and pages 182-183.
% For lndet information see: Ronald Barry and R. Kelley Pace, 
% "A Monte Carlo Estimator of the Log Determinant of Large Sparse Matrices", 
% Linear Algebra and its Applications", Volume 289, Number 1-3, 1999, pp. 41-54.
% and: R. Kelley Pace and Ronald P. Barry "Simulating Mixed Regressive
% Spatially autoregressive Estimators", 
% Computational Statistics, 1998, Vol. 13, pp. 397-418.
% ---------------------------------------------------

% written by:
% James P. LeSage, 1/2000
% Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial.econometrics.com

% NOTE: much of the speed for large problems comes from:
% the use of methods pioneered by Pace and Barry.
% R. Kelley Pace was kind enough to provide functions
% lndetmc, and lndetint from his spatial statistics toolbox
% for which I'm very grateful.


rflag = 0;
ldetflag = 1; % default to the fastest method
rflag = 0;
order = 50; liter = 30; % defaults
time1 = 0; 
time2 = 0;
time3 = 0;
results.order = order;
results.liter = liter;
hess_flag = 0;
lmin = -0.99;
lmax = 0.99;
rmin = -0.99;
rmax = 0.99;
    parm = [0.5
         0.5];

timet = clock; % start the clock for overall timing
% default options
options = optimset('fminsearch');

if nargin == 5
 if ~isstruct(info)
 error('sac: must supply the options as a structure variable');
 end;
options.MaxIter = 500;
 fields = fieldnames(info);
 nf = length(fields);
 for i=1:nf
    if strcmp(fields{i},'parm')
       parm = info.parm;
    elseif strcmp(fields{i},'convg')
       options.TolFun = info.convg;
    elseif strcmp(fields{i},'maxit')
        options.MaxIter  = info.maxit;
    elseif strcmp(fields{i},'rmin')
        rmin = info.rmin;
    elseif strcmp(fields{i},'rmax')
        rmax = info.rmax;
    elseif strcmp(fields{i},'lmin')
        lmin = info.lmin;
    elseif strcmp(fields{i},'lmax')
        lmax = info.lmax;
    elseif strcmp(fields{i},'hessian')
        hess_flag  = info.hessian;
    elseif strcmp(fields{i},'lflag')
        ldetflag = info.lflag;
    elseif strcmp(fields{i},'order')
        order = info.order;  results.order = order;
    elseif strcmp(fields{i},'iter')
    liter = info.iter; results.liter = liter;
    end;
 end;
elseif nargin == 4 % use default options
options = optimset('fminsearch');
else
 error('Wrong # of arguments to sac'); 
end; 


[n nvar] = size(x);
results.meth = 'sac';

[n1 n2] = size(W1);
if n1 ~= n2
error('sac: wrong size weight matrix W1');
elseif n1 ~= n
error('sac: wrong size weight matrix W1');
end;

[n1 n2] = size(W2);
if n1 ~= n2
error('sac: wrong size weight matrix W2');
elseif n1 ~= n
error('sac: wrong size weight matrix W2');
end;

results.y = y;      
results.nobs = n; results.nvar = nvar;
results.meth = 'sac';

% do lndet approximation calculations if needed
if ldetflag == 0 % no approximation
t0 = clock;    
out = lndetfull(W1,rmin,rmax);
time1 = etime(clock,t0);
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
det1 = [tt' outi];

t0 = clock;    
out = lndetfull(W2,lmin,lmax);
time1 = time1 + etime(clock,t0);
tt=lmin:.001:lmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
det2 = [tt' outi];

elseif ldetflag == 1 % use Pace and Barry, 1999 MC approximation

t0 = clock;    
out = lndetmc(order,liter,W1,rmin,rmax);
time1 = etime(clock,t0);
results.limit = [out.rho out.lo95 out.lndet out.up95];
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
det1 = [tt' outi];

t0 = clock;    
out = lndetmc(order,liter,W2,lmin,lmax);
time1 = time1 + etime(clock,t0);
results.limit = [out.rho out.lo95 out.lndet out.up95];
tt=lmin:.001:lmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
det2 = [tt' outi];

elseif ldetflag == 2 % use Pace and Barry, 1998 spline interpolation

t0 = clock;
out = lndetint(W1);
time1 = etime(clock,t0);
tt=.001:.001:1; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
det1 = [tt' outi];

t0 = clock;
out = lndetint(W2);
time1 = time1 + etime(clock,t0);
tt=.001:.001:1; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
det2 = [tt' outi];

end;

% find good starting values
% using Kelejian and Prucha GMM estimation
res0 = sac_gmm(y,x,W1,W2);
parm = [res0.rho
        res0.lam];

timeo = clock;
[pout,like,exitflag,output]=fminsearch('f_sac',parm,options,y,x,W1,W2,det1,det2);
time4 = etime(clock,timeo);

if exitflag == 0 
fprintf(1,'\n sac: convergence not obtained in %4d iterations \n',output.iterations);
end;
results.iter = output.iterations;

%results.lik = -like;
results.lik = -(n/2)-like;


rho = pout(1,1);
lam = pout(2,1);

% fill-in results
A = speye(n) - rho*sparse(W1);
B = speye(n) - lam*sparse(W2);
b0 = (B*x)\(B*A*y);
e = B*(A*y - x*b0);
results.beta = b0;
results.rho = rho;
results.lam = lam;
results.resid = e;
results.yhat = y-e;
sigu = e'*e;
sige = sigu/n;
results.sige = sige;

if (hess_flag == 0 & n <= 500)
% find asymptotic t-stats (from Anselin, 1982, pages 183-184
bhat = results.beta;
xpx = zeros(nvar+3,nvar+3);
BI = inv(B); AI = inv(A); WB = W2*BI; WA = W1*AI;
omeg = sige*eye(n); omegi = (1/sige)*eye(n);
% t-stats for beta
xpx(1:nvar,1:nvar) = (1/sige)*(x'*B'*B*x);
% t-stats for rho
%term1 = trace(WA.*WA);
term1 = trace(WA*WA);
term2 = trace(omeg*(B*WA*BI)'*omegi*(B*WA*BI));
term3 = (B*WA*x*bhat)'*omegi*(B*WA*x*bhat);
xpx(nvar+1,nvar+1) = term1+term2+term3;
% t-stats for lam
%term1 = trace(WB.*WB);
term1 = trace(WB*WB);
%term2 = (1/sige)*(W2*(A*y-x*bhat))'*(W2*(A*y-x*bhat));
term2 = trace(omeg*WB'*omegi*WB);
xpx(nvar+2,nvar+2) = term1+term2;
% sige,sige
xpx(nvar+3,nvar+3) = n/(2*sige*sige);
% off-diagonal terms bhat x rho
xpx(1:nvar,nvar+1) = (1/sige)*(x'*B'*B*WA*x*bhat);
xpx(nvar+1,1:nvar) = xpx(1:nvar,nvar+1)';
% off-diagonal terms bhat x lam
xpx(nvar+2,1:nvar) = zeros(1,nvar);
xpx(1:nvar,nvar+2) = xpx(nvar+2,1:nvar)';
% beta,sige = 0
% beta,rho = 0
% sige,rho
xpx(nvar+3,nvar+1) = (1/sige)*trace(W1*AI);
xpx(nvar+1,nvar+3) = xpx(nvar+3,nvar+1);
% sige,lambda
xpx(nvar+3,nvar+2) = (1/sige)*trace(W2*BI);
xpx(nvar+2,nvar+3) = xpx(nvar+3,nvar+2);
% off-diagonal terms rho x lam
% off-diagonal terms rho x lam
term1 = trace((WB)'*omegi*B*WA*BI*omeg);
term2 = trace(W2*WA*BI);
xpx(nvar+1,nvar+2) = term1+term2;
xpx(nvar+2,nvar+1) = xpx(nvar+1,nvar+2);
xpxi = invpd(xpx);
tmp = diag(abs(xpxi));
bvec = [results.beta
        results.rho
        results.lam];
results.tstat = bvec./sqrt(tmp(1:nvar+2,1));

elseif (n > 500 | hess_flag == 1) % use numerical hessian
t0 = clock;

parm = [results.beta
        results.rho
        results.lam
        results.sige];

if ldetflag == 0
hessn = hessian('f2_sac',parm,y,x,W1,W2,det1,det2);
else
hessn = hessian('f2_sac',parm,y,x,W1,W2,det1,det2);
end;

xpxi = invpd(-hessn);
xpxi = abs(diag(xpxi(1:nvar+2,1:nvar+2)));
tmp = [results.beta
       results.rho
       results.lam];
results.tstat = tmp./sqrt(xpxi);
time3 = etime(clock,t0);

end;    
    
% r-squared and rbar-squared
ym = y - mean(y);
rsqr1 = sigu;
rsqr2 = ym'*ym;
results.rsqr = 1.0 - rsqr1/rsqr2; % r-squared
rsqr1 = rsqr1/(n-nvar);
rsqr2 = rsqr2/(n-1.0);
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
results.time = etime(clock,timet);
results.time1 = time1;
results.time2 = time2;
results.time3 = time3;
results.time4 = time4;


