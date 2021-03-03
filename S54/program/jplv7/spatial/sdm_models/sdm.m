function results = sdm(y,x,W,info)
% PURPOSE: computes spatial durbin model estimates
%         (I-rho*W)y = a + X*B1 + W*X*B2 + e, using sparse algorithms
% ---------------------------------------------------
%  USAGE: results = sdm(y,x,W,info)
%  where: y = dependent variable vector
%         x = independent variables matrix
%         W = contiguity matrix (standardized)
%       info = an (optional) structure variable with input options:
%       info.rmin = (optional) minimum value of rho to use in search  
%       info.rmax = (optional) maximum value of rho to use in search    
%       info.convg = (optional) convergence criterion (default = 1e-8)
%       info.maxit = (optional) maximum # of iterations (default = 500)
%       info.lflag = 0 for full lndet computation (default = 1, fastest)
%                  = 1 for MC lndet approximation (fast for very large problems)
%                  = 2 for Spline lndet approximation (medium speed)
%       info.order = order to use with info.lflag = 1 option (default = 50)
%       info.iter  = iterations to use with info.lflag = 1 option (default = 30) 
%       info.lndet = a matrix returned by sar, sar_g, sarp_g, etc.
%                    containing log-determinant information to save time
% ---------------------------------------------------
%  RETURNS: a structure
%         results.meth  = 'sdm'
%         results.beta  = bhat [a B1 B2]' a k+(k-1) x 1 vector
%         results.rho   = rho 
%         results.tstat = t-statistics (last entry is rho)
%         results.yhat  = yhat
%         results.resid = residuals
%         results.sige  = sige
%         results.rsqr  = rsquared
%         results.rbar  = rbar-squared
%         results.lik   = log likelihood
%         results.nobs  = nobs
%         results.nvar  = nvars (includes lam)
%         results.y     = y data vector
%         results.iter  = # of iterations taken
%         results.rmax  = 1/max eigenvalue of W (or rmax if input)
%         results.rmin  = 1/min eigenvalue of W (or rmin if input)
%         results.lflag = lflag from input
%         results.miter = info.iter option from input
%         results.order = info.order option from input
%         results.limit = matrix of [rho lower95,logdet approx, upper95] intervals
%                         for the case of lflag = 1
%         results.time1 = time for log determinant calcluation
%         results.time2 = time for eigenvalue calculation
%         results.time3 = time for hessian or information matrix calculation
%         results.time4 = time for optimization
%         results.time  = total time taken       
%         results.lndet = a matrix containing log-determinant information
%                          (for use in later function calls to save time)
%  --------------------------------------------------
%  SEE ALSO: prt(results)
% ---------------------------------------------------
%  NOTES: constant term should be in 1st column of the x-matrix
%         constant is excluded from B2 estimates
%        YOU MUST have a constant term in the model for this code to work
%  if you use lflag = 1 or 2, info.rmin will be set = -1 
%                             info.rmax will be set = 1
%  For n < 1000 you should use lflag = 0 to get exact results  
% ---------------------------------------------------
% REFERENCES: R. Kelley Pace and R. Barry, 
% ``Simulating mixed regressive spatial autoegressive estimators'', 
% Computational Statistics, 1998 Vol. 13, pp. 397-418.
% For lndet information see: Ronald Barry and R. Kelley Pace, 
% "A Monte Carlo Estimator of the Log Determinant of Large Sparse Matrices", 
% Linear Algebra and its Applications", Volume 289, Number 1-3, 1999, pp. 41-54.
% ---------------------------------------------------

% written by:
% James P. LeSage, 4/2002, last updated 7/2003
% Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

% NOTE: much of the speed for large problems comes from:
% the use of methods pioneered by Pace and Barry.
% R. Kelley Pace was kind enough to provide functions
% lndetmc, and lndetint from his spatial statistics toolbox
% for which I'm very grateful.


% check if the user handled the intercept term okay
n = length(y);
if sum(x(:,1)) ~= n
tst = sum(x); % we may have no intercept term
ind = find(tst == n); % we do have an intercept term
 if length(ind) > 0
 error('sdm: intercept term must be in first column of the x-matrix');
 elseif length(ind) == 0 % case of no intercept term
 xsdm = [x W*x];
 cflag = 0;
 end;
elseif sum(x(:,1)) == n % we have an intercept in the right place
 xsdm = [x W*x(:,2:end)];
 cflag = 1;
end;

[nobs,nvar] = size(x);

% just call sar function
if nargin == 4
results = sar(y,xsdm,W,info);
elseif nargin == 3
results = sar(y,xsdm,W);
else
error('sdm: wrong # of input arguments to sdm');
end;

results.meth = 'sdm';
results.nvar = nvar;
results.cflag = cflag;
