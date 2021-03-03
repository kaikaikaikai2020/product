function results=sar_gmm(y,x,W)
% PURPOSE: computes Generalized Moments Estimates for Spatial Autoregressive Model
%           y = rho*W*y + XB + e
% ---------------------------------------------------
%  USAGE: results = sar_gmm(y,x,W)
%  where: y = dependent variable vector
%         x = independent variables matrix
%             (with intercept vector in the 1st column of x)
%         W = sparse contiguity matrix (standardized)
% ---------------------------------------------------
%  RETURNS: a structure
%         results.meth       = 'sar_gmm'
%         results.beta       = bhat
%         results.tstat      = asymp t-stats
%         results.rho        = rho
%         results.rhotstat   = t-stat of rho (under normality assumption)
%         results.GMsige     = GM-estimated variance
%         results.yhat       = yhat = B*x*bhat, B=inv(I - rho*W)
%         results.resid      = residuals, y - yhat
%         results.sige       = sige = e'*e/n, e = y - yhat
%         results.rsqr       = rsquared
%         results.rbar       = rbar-squared
%         results.se         = Standard errors from EGLS
%         results.nobs       = number of observations
%         results.nvar       = number of variables 
%         results.time       = total time taken
% ---------------------------------------------------
% %  SEE ALSO: prt(results), sar, sar_g
% ---------------------------------------------------
% REFERENCES: Luc Anselin Spatial Econometrics (1988) pages 182-183.
% Kelejian, H., and  Prucha, I.R.  (1998). A Generalized Spatial Two-Stage
% Least Squares Procedure for Estimating a Spatial Autoregressive
% Model with Autoregressive Disturbances. Journal of Real
% Estate and Finance Economics,  17, 99-121.
% ---------------------------------------------------

% written by: Jim LeSage
% Adapted from Shawn Bucholtz code for the SEM model case

timet = clock; % start the clock for overall timing

% error checking on inputs
xsum = sum(x);
[n,k] = size(x);
ind = find(xsum == n);
iflag = 0;
if length(ind) > 0 % we have an intercept
    if ind ~= 1
    warning('intercept must be in 1st column of the x-matrix');
    end;
    iflag = 1;
end;

results.meth = 'sar_gmm';
time1 = 0; 
time2 = 0;



%Estimated OLS to get a vector of residuals
[n, nvar]=size(x);
results.nobs=n;
results.nvar=nvar;
    if iflag == 1
       Wy = sparse(W)*y;
       Wx = sparse(W)*x(:,2:end);
        z = [x Wx W*Wx];
       o1 = tsls(y,Wy,x,z);
      rho = o1.beta(1,1);
     bhat = o1.beta(2:end,1);
 rhotstat = o1.tstat(1,1);
   btstat = o1.tstat(2:end,1);
     sige = (o1.resid'*o1.resid)/n;
    elseif iflag == 0
    Wy = sparse(W)*y;
    Wx = sparse(W)*x;
    z = [x Wx W*Wx];
    o1 = tsls(y,Wy,x,z);
    rho = o1.beta(1,1);
    bhat = o1.beta(2:end,1);
    rhotstat = o1.tstat(1,1);
    btstat = o1.tstat(2:end,1);
    sige = (o1.resid'*o1.resid)/n;
    end;    
    
results.rho = rho;
results.beta = bhat;
results.sige = sige;
results.tstat = btstat;
results.rhotstat = rhotstat;

sigu = results.sige*n;
ym = y - mean(y);
rsqr1 = sigu;
rsqr2 = ym'*ym;
results.rsqr = 1.0 - rsqr1/rsqr2; % r-squared
rsqr1 = rsqr1/(n-nvar);
rsqr2 = rsqr2/(n-1.0);
if rsqr2 ~= 0
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
else
    results.rbar = results.rsqr;
end;

time2 = etime(clock,timet);

results.time = time2;
