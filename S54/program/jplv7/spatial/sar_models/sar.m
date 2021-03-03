function results = sar(y,x,W,info)
% PURPOSE: computes spatial autoregressive model estimates
%           y = p*W*y + X*b + e, using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE: results = sar(y,x,W,info)
%  where:  y = dependent variable vector
%          x = explanatory variables matrix
%          W = standardized contiguity matrix 
%       info = an (optional) structure variable with input options:
%       info.rmin  = (optional) minimum value of rho to use in search (default = -1) 
%       info.rmax  = (optional) maximum value of rho to use in search (default = +1)    
%       info.eig   = 0 for default rmin = -1,rmax = +1, 1 for eigenvalue calculation of these
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
%         results.meth  = 'sar'
%         results.beta  = bhat (nvar x 1) vector
%         results.rho   = rho
%         results.tstat = asymp t-stat (last entry is rho)
%         results.bstd  = std of betas (nvar x 1) vector
%         results.pstd  = std of rho
%         results.yhat  = yhat         (nobs x 1) vector
%         results.resid = residuals    (nobs x 1) vector
%         results.sige  = sige = (y-p*W*y-x*b)'*(y-p*W*y-x*b)/n
%         results.rsqr  = rsquared
%         results.rbar  = rbar-squared
%         results.lik   = log likelihood
%         results.nobs  = # of observations
%         results.nvar  = # of explanatory variables in x 
%         results.y     = y data vector
%         results.iter  = # of iterations taken
%         results.rmax  = 1/max eigenvalue of W (or rmax if input)
%         results.rmin  = 1/min eigenvalue of W (or rmin if input)
%         results.lflag = lflag from input
%         results.liter = info.iter option from input
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
% --------------------------------------------------
%  NOTES: if you use lflag = 1 or 2, info.rmin will be set = -1 
%                                    info.rmax will be set = 1
%         For n < 1000 you should use lflag = 0 to get exact results                                    
% --------------------------------------------------  
%  SEE ALSO: prt(results), sac, sem, sdm, sar, far
% ---------------------------------------------------
% REFERENCES: Anselin (1988), pages 180-182.
% For lndet information see: Ronald Barry and R. Kelley Pace, 
% "A Monte Carlo Estimator of the Log Determinant of Large Sparse Matrices", 
% Linear Algebra and its Applications", Volume 289, Number 1-3, 1999, pp. 41-54.
% and: R. Kelley Pace and Ronald P. Barry 
% "Simulating Mixed Regressive Spatially autoregressive Estimators", 
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


time1 = 0; 
time2 = 0;
time3 = 0;
time4 = 0;

timet = clock; % start the clock for overall timing

% check size of user inputs for comformability
[n nvar] = size(x); [n1 n2] = size(W);
if n1 ~= n2
error('sar: wrong size weight matrix W');
elseif n1 ~= n
error('sar: wrong size weight matrix W');
end;
[nchk junk] = size(y);
if nchk ~= n
error('sar: wrong size y vector input');
end;

% if we have no options, invoke defaults
if nargin == 3
    info.lflag = 1;
end;


% parse input options
[rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,miter,options] = sar_parse(info);

% compute eigenvalues or limits
[rmin,rmax,time2] = sar_eigs(eflag,W,rmin,rmax,n);

% do log-det calculations
[detval,time1] = sar_lndet(ldetflag,W,rmin,rmax,detval,order,miter);

t0 = clock;
          Wy = sparse(W)*y;
          AI = x'*x;
          b0 = AI\(x'*y);
          bd = AI\(x'*Wy);
          e0 = y - x*b0;
          ed = Wy - x*bd;
          epe0 = e0'*e0;
          eped = ed'*ed;
          epe0d = ed'*e0;

% step 1) do regressions
% step 2) maximize concentrated likelihood function;
    options = optimset('fminbnd');
    [p,liktmp,exitflag,output] = fminbnd('f_sar',rmin,rmax,options,detval,epe0,eped,epe0d,n);
   
time4 = etime(clock,t0);

if exitflag == 0
fprintf(1,'\n sar: convergence not obtained in %4d iterations \n',output.iterations);
end;
results.iter = output.iterations;

% step 3) find b,sige maximum likelihood estimates
results.beta = b0 - p*bd; 
results.rho = p; 
bhat = results.beta;
results.sige = (1/n)*(e0-p*ed)'*(e0-p*ed); 
sige = results.sige;

e = (e0 - p*ed);
yhat = (speye(n) - p*W)\(x*bhat);
results.yhat = yhat;
results.resid = y - yhat;

parm = [results.beta
        results.rho
        results.sige];

results.lik = f2_sar(parm,y,x,W,detval);

if n <= 500
t0 = clock;
% asymptotic t-stats based on information matrix
% (page 80-81 Anselin, 1980)
B = eye(n) - p*W; 
BI = inv(B); WB = W*BI;
pterm = trace(WB*WB + WB*WB');
xpx = zeros(nvar+2,nvar+2);               % bhat,bhat
xpx(1:nvar,1:nvar) = (1/sige)*(x'*x);     % bhat,rho
xpx(1:nvar,nvar+1) = (1/sige)*x'*W*BI*x*bhat;
xpx(nvar+1,1:nvar) = xpx(1:nvar,nvar+1)'; % rho,rho
xpx(nvar+1,nvar+1) = (1/sige)*bhat'*x'*BI'*W'*W*BI*x*bhat + pterm;
xpx(nvar+2,nvar+2) = n/(2*sige*sige);     %sige,sige
xpx(nvar+1,nvar+2) = (1/sige)*trace(WB);  % rho,sige
xpx(nvar+2,nvar+1) = xpx(nvar+1,nvar+2);
xpxi = invpd(xpx);

tmp = diag(abs(xpxi(1:nvar+1,1:nvar+1)));
bvec = [results.beta
        results.rho];
tmps = bvec./(sqrt(tmp));
results.tstat = tmps;
results.bstd = sqrt(tmp(1:nvar,1));
results.pstd = sqrt(tmp(nvar,1));
time3 = etime(clock,t0);

else  % asymptotic t-stats using numerical hessian
t0 = clock;
% just computes the diagonal
dhessn = hessian('f2_sar',parm,y,x,W,detval);
hessi = invpd(dhessn);

tvar = abs(diag(hessi));
tmp = [results.beta
       results.rho];
results.tstat = tmp./sqrt(tvar(1:end-1,1));
results.bstd = sqrt(tvar(1:end-2,1));
results.pstd = sqrt(tvar(end-1,1));
time3 = etime(clock,t0);

end; % end of t-stat calculations


ym = y - mean(y);       % r-squared, rbar-squared
rsqr1 = results.resid'*results.resid;
rsqr2 = ym'*ym;
results.rsqr = 1.0-rsqr1/rsqr2;   % r-squared
rsqr1 = rsqr1/(n-nvar);
rsqr2 = rsqr2/(n-1.0);

% return stuff
results.meth = 'sar';
results.y = y;      
results.nobs = n; 
results.nvar = nvar;
results.rmax = rmax;      
results.rmin = rmin;
results.lflag = ldetflag;
results.order = order;
results.miter = miter;
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
results.time = etime(clock,timet);
results.time1 = time1;
results.time2 = time2;
results.time3 = time3;
results.time4 = time4;
results.lndet = detval;



function llike = f_sar(rho,detval,epe0,eped,epe0d,n)
% PURPOSE: evaluates concentrated log-likelihood for the 
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:llike = f_sar(rho,detval,epe0,eped,epe0d,n)
%  where: rho  = spatial autoregressive parameter
%         detval = a matrix with vectorized log-determinant information
%         epe0   = see below
%         eped   = see below
%         eoe0d  = see below
%         n      = # of obs
%          b0 = AI*xs'*ys;
%          bd = AI*xs'*Wys;
%          e0 = ys - xs*b0;
%          ed = Wys - xs*bd;
%          epe0 = e0'*e0;
%          eped = ed'*ed;
%          epe0d = ed'*e0;
% ---------------------------------------------------
%  RETURNS: a  scalar equal to minus the log-likelihood
%           function value at the parameter rho
% ---------------------------------------------------                         
%  NOTE: this is really two functions depending
%        on nargin = 3 or nargin = 4 (see the function)
%  --------------------------------------------------
%  SEE ALSO: sar, f_far, f_sac, f_sem
% ---------------------------------------------------

% written by: James P. LeSage 1/2000
% University of Toledo
% Department of Economics
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

if nargin == 6
gsize = detval(2,1) - detval(1,1);
% Note these are actually log detvalues
i1 = find(detval(:,1) <= rho + gsize);
i2 = find(detval(:,1) <= rho - gsize);
i1 = max(i1);
i2 = max(i2);
index = round((i1+i2)/2);
if isempty(index)
index = 1;
end;

detm = detval(index,2); 

z = epe0 - 2*rho*epe0d + rho*rho*eped;

llike = (n/2)*log(z) - detm;

else

error('f_sar: Wrong # of input arguments');

end;

function llike = f2_sar(parm,y,x,W,detval)
% PURPOSE: evaluates log-likelihood -- given ML estimates
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:llike = f2_sar(parm,y,X,W,ldet)
%  where: parm = vector of maximum likelihood parameters
%                parm(1:k-2,1) = b, parm(k-1,1) = rho, parm(k,1) = sige
%         y    = dependent variable vector (n x 1)
%         X    = explanatory variables matrix (n x k)
%         W    = spatial weight matrix
%         ldet = matrix with [rho log determinant] values
%                computed in sar.m using one of Kelley Pace's routines  
% ---------------------------------------------------
%  RETURNS: a  scalar equal to minus the log-likelihood
%           function value at the ML parameters
%  --------------------------------------------------
%  NOTE: this is really two functions depending
%        on nargin = 4 or nargin = 5 (see the function)
% ---------------------------------------------------
%  SEE ALSO: sar, f2_far, f2_sac, f2_sem
% ---------------------------------------------------

% written by: James P. LeSage 1/2000
% University of Toledo
% Department of Economics
% Toledo, OH 43606
% jlesage@spatial.econometrics.com

n = length(y); 
k = length(parm);
b = parm(1:k-2,1);
rho = parm(k-1,1);
sige = parm(k,1);

gsize = detval(2,1) - detval(1,1);
i1 = find(detval(:,1) <= rho + gsize);
i2 = find(detval(:,1) <= rho - gsize);
i1 = max(i1);
i2 = max(i2);
index = round((i1+i2)/2);
if isempty(index)
index = 1;
end;
detm = detval(index,2);

e = y-x*b-rho*sparse(W)*y;
epe = e'*e;
tmp2 = 1/(2*sige);
llike = -(n/2)*log(pi) - (n/2)*log(sige) + detm - tmp2*epe;


function [rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,iter,options] = sar_parse(info)
% PURPOSE: parses input arguments for sar model
% ---------------------------------------------------
%  USAGE: [rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,iter,options] = sar_parse(info)
% where info contains the structure variable with inputs 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------

% set defaults
options = zeros(1,18); % optimization options for fminbnd
options(1) = 0; 
options(2) = 1.e-6; 
options(14) = 500;

eflag = 0;     % default to not computing eigenvalues
ldetflag = 1;  % default to 1999 Pace and Barry MC determinant approx
order = 50;    % there are parameters used by the MC det approx
iter = 30;     % defaults based on Pace and Barry recommendation
rmin = -1;     % use -1,1 rho interval as default
rmax = 1;
detval = 0;    % just a flag
convg = 0.0001;
maxit = 500;

fields = fieldnames(info);
nf = length(fields);
if nf > 0
    
 for i=1:nf
    if strcmp(fields{i},'rmin')
        rmin = info.rmin;  eflag = 0;
    elseif strcmp(fields{i},'rmax')
        rmax = info.rmax; eflag = 0;
    elseif strcmp(fields{i},'convg')
        options(2) = info.convg;
    elseif strcmp(fields{i},'maxit')
        options(14) = info.maxit;  
    elseif strcmp(fields{i},'lndet')
    detval = info.lndet;
    ldetflag = -1;
    eflag = 0;
    rmin = detval(1,1);
    nr = length(detval);
    rmax = detval(nr,1);
    elseif strcmp(fields{i},'lflag')
        tst = info.lflag;
        if tst == 0,
        ldetflag = 0; % compute full lndet, no approximation
        elseif tst == 1,
        ldetflag = 1; % use Pace-Barry approximation
        elseif tst == 2,
        ldetflag = 2; % use spline interpolation approximation
        else
        error('sar: unrecognizable lflag value on input');
        end;
    elseif strcmp(fields{i},'order')
        order = info.order;  
    elseif strcmp(fields{i},'eig')
        eflag = info.eig;  
    elseif strcmp(fields{i},'iter')
        iter = info.iter; 
    end;
 end;
 
else, % the user has input a blank info structure
      % so we use the defaults
end; 

function [rmin,rmax,time2] = sar_eigs(eflag,W,rmin,rmax,n);
% PURPOSE: compute the eigenvalues for the weight matrix
% ---------------------------------------------------
%  USAGE: [rmin,rmax,time2] = far_eigs(eflag,W,rmin,rmax,W)
% where eflag is an input flag, W is the weight matrix
%       rmin,rmax may be used as default outputs
% and the outputs are either user-inputs or default values
% ---------------------------------------------------


if eflag == 1 % do eigenvalue calculations
t0 = clock;
opt.tol = 1e-3; opt.disp = 0;
lambda = eigs(sparse(W),speye(n),1,'SR',opt);  
rmin = real(1/lambda);   
rmax = 1.0;
time2 = etime(clock,t0);
else % use rmin,rmax arguments from input or defaults -1,1
time2 = 0;
end;


function [detval,time1] = sar_lndet(ldetflag,W,rmin,rmax,detval,order,iter);
% PURPOSE: compute the log determinant |I_n - rho*W|
% using the user-selected (or default) method
% ---------------------------------------------------
%  USAGE: detval = far_lndet(lflag,W,rmin,rmax)
% where eflag,rmin,rmax,W contains input flags 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------


% do lndet approximation calculations if needed
if ldetflag == 0 % no approximation
t0 = clock;    
out = lndetfull(W,rmin,rmax);
time1 = etime(clock,t0);
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
detval = [tt' outi];
    
elseif ldetflag == 1 % use Pace and Barry, 1999 MC approximation

t0 = clock;    
out = lndetmc(order,iter,W,rmin,rmax);
time1 = etime(clock,t0);
results.limit = [out.rho out.lo95 out.lndet out.up95];
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
detval = [tt' outi];

elseif ldetflag == 2 % use Pace and Barry, 1998 spline interpolation

t0 = clock;
out = lndetint(W,rmin,rmax);
time1 = etime(clock,t0);
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
detval = [tt' outi];

elseif ldetflag == -1 % the user fed down a detval matrix
    time1 = 0;
        % check to see if this is right
        if detval == 0
            error('sar: wrong lndet input argument');
        end;
        [n1,n2] = size(detval);
        if n2 ~= 2
            error('sar: wrong sized lndet input argument');
        elseif n1 == 1
            error('sar: wrong sized lndet input argument');
        end;          
end;
