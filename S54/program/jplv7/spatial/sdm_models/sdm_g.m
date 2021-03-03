function results = sdm_g(y,x,W,ndraw,nomit,prior)
% PURPOSE: Bayesian estimates of the heteroscedastic spatial durbin model
%         (I-rho*W)y = a + X*B1 + W*X*B2 + e, e = N(0,sige*V), V = diag(v1,v2,...vn)
%          r/vi = ID chi(r)/r, r = Gamma(m,k)
%          a, B1, B2 = diffuse
%          1/sige = Gamma(nu,d0), 
%          rho = Uniform(rmin,rmax), or rho = beta(a1,a2); 
%-------------------------------------------------------------
% USAGE: results = sdm_g(y,x,W,ndraw,nomit,prior)
% where: y = dependent variable vector (nobs x 1)
%        x = independent variables matrix (nobs x nvar), with constant term in 1st column
%        W = spatial weight matrix (standardized, row-sums = 1)
%    ndraw = # of draws
%    nomit = # of initial draws omitted for burn-in            
%    prior = a structure variable with:
%            prior.rval  = r prior hyperparameter, default=4
%            prior.novi  = 1 turns off sampling for vi, producing homoscedastic model            
%            prior.m     = informative Gamma(m,k) prior on r
%            prior.k     = (default: not used)
%            prior.nu    = informative Gamma(nu,d0) prior on sige
%            prior.d0    = default: nu=0,d0=0 (diffuse prior)
%            prior.rmin  = (optional) min rho used in sampling (default = 0)
%            prior.rmax  = (optional) max rho used in sampling (default = 1)  
%            prior.a1    = parameter for beta(a1,a2) prior on rho (default = 1.01)
%            prior.a2    = (default = 1.01) see: 'help beta_prior'
%            prior.lflag = 0 for full lndet computation (default = 1, fastest)
%                        = 1 for MC approx (fast for large problems)
%                        = 2 for Spline approx (medium speed) 
%            prior.order = order to use with prior.lflag = 1 option (default = 50)
%            prior.iter  = iters to use with prior.lflag = 1 option (default = 30)   
%            prior.lndet = a matrix returned by sar, sar_g, sarp_g, etc.
%                          containing log-determinant information to save time
%            prior.dflag = 0 for numerical integration, 1 for Metropolis-Hastings
%                          (default = 0)
%-------------------------------------------------------------
% RETURNS:  a structure:
%          results.meth   = 'sdm_g'
%          results.beta   = posterior mean of bhat (based on draws)
%          results.rho    = posterior mean of rho  (based on draws)
%          results.sige   = posterior mean of sige (based on draws)
%          results.sigma  = sige based on (e'*e)/(n-k)
%          results.bdraw  = bhat draws (ndraw-nomit x nvar)
%          results.pdraw  = rho  draws (ndraw-nomit x 1)
%          results.sdraw  = sige draws (ndraw-nomit x 1)
%          results.vmean  = mean of vi draws (nobs x 1) 
%          results.rdraw  = r draws (ndraw-nomit x 1) (if m,k input)
%          results.mlike  = log marginal posterior, a scalar
%          results.r      = value of hyperparameter r (if input)
%          results.nobs   = # of observations
%          results.nvar   = # of variables in x-matrix
%          results.ndraw  = # of draws
%          results.nomit  = # of initial draws omitted
%          results.y      = y-vector from input (nobs x 1)
%          results.yhat   = mean of posterior predicted (nobs x 1)
%          results.resid  = residuals, based on posterior means
%          results.rsqr   = r-squared based on posterior means
%          results.rbar   = adjusted r-squared
%          results.nu     = nu prior parameter
%          results.d0     = d0 prior parameter
%          results.a1     = a1 parameter for beta prior on rho from input, or default value
%          results.a2     = a2 parameter for beta prior on rho from input, or default value
%          results.time1  = time for eigenvalue calculation
%          results.time2  = time for log determinant calcluation
%          results.time3  = time for sampling
%          results.time   = total time taken  
%          results.rmax   = 1/max eigenvalue of W (or rmax if input)
%          results.rmin   = 1/min eigenvalue of W (or rmin if input)          
%          results.tflag  = 'plevel' (default) for printing p-levels
%                         = 'tstat' for printing bogus t-statistics 
%          results.lflag  = lflag from input
%          results.dflag  = the value of dflag from input
%          results.iter   = prior.iter option from input
%          results.order  = prior.order option from input
%          results.limit  = matrix of [rho lower95,logdet approx, upper95] 
%                           intervals for the case of lflag = 1
%          results.lndet = a matrix containing log-determinant information
%                          (for use in later function calls to save time)
%          results.novi  = novi from input (or default)
%          results.acc   = acceptance rate for M-H sampling
%          results.mlike = log marginal likelihood (a vector ranging over
%                          rho values that can be integrated for model comparison)
% --------------------------------------------------------------
% NOTES: constant term should be in 1st column of the x-matrix
%        constant is excluded from B2 estimates
%        YOU MUST have a constant term in the model for this code to work
% - use either improper prior.rval 
%          or informative Gamma prior.m, prior.k, not both of them
% - for n < 1000 you should use lflag = 0 to get exact results  
% - If you wish to place informative priors on beta, use sar_g function
%   and add W*x as explanatory variables with diffuse or informative priors
% - use a1 = 1.0 and a2 = 1.0 for uniform prior on rho
% --------------------------------------------------------------
% SEE ALSO: (sdm_gd, sdm_gd2 demos) sdm_gc.m, prt
% --------------------------------------------------------------
% REFERENCES: James P. LeSage, `Bayesian Estimation of Spatial Autoregressive
%             Models',  International Regional Science Review, 1997 
%             Volume 20, number 1\&2, pp. 113-129.
% also, R. Kelley Pace and Ronald P. Barry 
% "Simulating Mixed Regressive Spatially autoregressive Estimators", 
% Computational Statistics, 1998, Vol. 13, pp. 397-418.
%----------------------------------------------------------------

% written by:
% James P. LeSage, 12/2001, updated 7/2003
% Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com


% NOTE: some of the speed for large problems comes from:
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

if nargin == 6
results = sar_g(y,xsdm,W,ndraw,nomit,prior);
elseif nargin == 5
results = sar_g(y,xsdm,W,ndraw,nomit);
else
error('sdm_g: wrong # of input arguments');
end;

results.meth = 'sdm_g';
results.nvar = nvar;
results.cflag = cflag;
