function results = sem_panel(y,x,W,T,info)
% PURPOSE: computes spatial error model estimates for spatial panels (N regions*T time periods)
%           y = XB + u,  u = p*W*u + e, using sparse algorithms
% Supply data sorted first by time and then by spatial units, so first region 1,
% region 2, et cetera, in the first year, then region 1, region 2, et
% cetera in the second year, and so on
% sem_panel computes y and x in deviation of the spatial and/or time means
% (see Baltagi, 2001, Econometric Analysis of Panel Data, ch. 2 and ch. 3)
% ---------------------------------------------------
%  USAGE: results = sem_panel(y,x,W,T,info)
%  where: y = dependent variable vector
%         x = independent variables matrix 
%         W = spatial weights matrix (standardized)
%         T = number of points in time
%       info = an (optional) structure variable with input options:
%       info.model = 0 pooled model without fixed effects (default, x may contain an intercept)
%                  = 1 spatial fixed effects (x may not contain an intercept)
%                  = 2 time period fixed effects (x may not contain an intercept)
%                  = 3 spatial and time period fixed effects (x may not contain an intercept)
%       info.rmin  = (optional) minimum value of rho to use in search  
%       info.rmax  = (optional) maximum value of rho to use in search    
%       info.convg = (optional) convergence criterion (default = 1e-4)
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
%         results.meth  = 'psem' if infomodel=0
%                       = 'semsfe' if info.model=1
%                       = 'semtfe' if info.model=2
%                       = 'semstfe' if info.model=3
%         results.beta  = bhat
%         results.rho   = rho (p above)
%         results.tstat = asymp t-stats (last entry is rho=spatial autocorrelation coefficient)
%         results.yhat  = yhat
%         results.resid = residuals
%         results.sige  = sige = e'(I-p*W)'*(I-p*W)*e/nobs
%         results.rsqr  = rsquared
%         results.rbar  = rbarsquared
%         results.sfe   = spatial fixed effects (if info.model=1 or 3)
%         results.tfe   = time period fixed effects (if info.model=2 or 3)
%         results.con   = intercept (if info.model=3)
%         results.lik   = log likelihood
%         results.nobs  = # of observations
%         results.nvar  = # of explanatory variables in x
%         results.tnvar = nvar + # fixed effects
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
%  --------------------------------------------------
%  NOTES: if you use lflag = 1 or 2, info.rmin will be set = -1 
%                                    info.rmax will be set = 1
%         For number of spatial units < 500 you should use lflag = 0 to get exact results
% ---------------------------------------------------
%
% written by: J.Paul Elhorst 11/2004
% University of Groningen
% Department of Economics
% 9700AV Groningen
% the Netherlands
% j.p.elhorst@eco.rug.nl
%
% REFERENCES: 
% "Specification and Estimation of Spatial Panel Data Models",
% International Regional Science Review, Vol. 26, pp. 244-268.
% Formulas for information matrix are not in this paper, I derived them
% later

% This function is based on James. P LeSage's function SEM

time1 = 0; 
time2 = 0;
time3 = 0;

timet = clock; % start the clock for overall timing

% if we have no options, invoke defaults
if nargin == 4
    info.lflag = 1;
    info.model = 0;
    fprintf(1,'default: pooled model without fixed effects \n');
end;

fields = fieldnames(info);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'model') model = info.model;
        end
    end
end
if model==0
    results.meth='psem';
elseif model==1
    results.meth='semsfe';
elseif model==2
    results.meth='semtfe';
elseif model==3
    results.meth='semstfe';
else
    error('sem_panel: wrong input number of info.model');
end

% check size of user inputs for comformability
[nobs nvar] = size(x);
[N Ncol] = size(W);
if N ~= Ncol
error('sem: wrong size weight matrix W');
elseif N ~= nobs/T
error('sem: wrong size weight matrix W or matrix x');
end;
[nchk junk] = size(y);
if nchk ~= nobs
error('sar: wrong size vector y or matrix x');
end;

% return the easy stuff
results.y = y;
results.nobs = nobs;
results.nvar = nvar; 

% parse input options
[rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,miter,options] = sem_parse(info); %function of LeSage

% compute eigenvalues or limits
[rmin,rmax,time2] = sem_eigs(eflag,W,rmin,rmax,N); %function of LeSage

results.rmin = rmin;
results.rmax = rmax;
results.lflag = ldetflag;
results.miter = miter;
results.order = order;

% do log-det calculations
[detval,time1] = sem_lndet(ldetflag,W,rmin,rmax,detval,order,miter); % function of LeSage

% demeaning of the y and x variables, depending on (info.)model

if (model==1 | model==3);
meanny=zeros(N,1);
meannx=zeros(N,nvar);
for i=1:N
    ym=zeros(T,1);
    xm=zeros(T,nvar);
    for t=1:T
        ym(t)=y(i+(t-1)*N,1);
        xm(t,:)=x(i+(t-1)*N,:);
    end
    meanny(i)=mean(ym);
    meannx(i,:)=mean(xm);
end
clear ym xm;
end % if statement

if ( model==2 | model==3)
meanty=zeros(T,1);
meantx=zeros(T,nvar);
for i=1:T
    t1=1+(i-1)*N;t2=i*N;
    ym=y([t1:t2],1);
    xm=x([t1:t2],:);
    meanty(i)=mean(ym);
    meantx(i,:)=mean(xm);
end
clear ym xm;
end % if statement
    
en=ones(T,1);
et=ones(N,1);
ent=ones(nobs,1);

if model==1;
    ywith=y-kron(en,meanny);
    xwith=x-kron(en,meannx);
elseif model==2
    ywith=y-kron(et,meanty);
    xwith=x-kron(et,meantx);
elseif model==3
    ywith=y-kron(en,meanny)-kron(et,meanty)+kron(ent,mean(y));
    xwith=x-kron(en,meannx)-kron(et,meantx)+kron(ent,mean(x));
else
    ywith=y; 
    xwith=x;
end % if statement
    
t0 = clock;

for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    Wx([t1:t2],:)= sparse(W)*xwith([t1:t2],:);
    Wy([t1:t2],1)= sparse(W)*ywith([t1:t2],1);
end

options = optimset('MaxIter',maxit);
rho = 0.5;
converge = 1;
criteria = 1e-4;
iter = 1;

while (converge > criteria) & (iter < maxit)

xs = xwith - rho*Wx;
ys = ywith - rho*Wy;
b = (xs'*xs)\(xs'*ys);
e = (ywith - xwith*b);
rold = rho;
[rho,like,exitflag,output] = fminbnd('f_sempanel',rmin,rmax,options,e,W,detval,T);
converge = abs(rold - rho);
iter = iter + 1;
end;


liktmp = like;
time4 = etime(clock,t0);

if exitflag == maxit
fprintf(1,'\n sem: convergence not obtained in %4d iterations \n',output.iterations);
end;
results.iter = output.iterations;

% compute results 
results.beta= b;
results.rho = rho;

if model==1
    results.sfe=meanny-meannx*results.beta; % including intercept
    results.yhat=x*results.beta+kron(en,results.sfe);
    tnvar=nvar+N;
elseif model==2
    results.tfe=meanty-meantx*results.beta; % including intercept
    results.yhat=x*results.beta+kron(et,results.tfe);
    tnvar=nvar+T;
elseif model==3
    intercept=mean(y)-mean(x)*results.beta; % intercept calculated separately
    results.con=intercept;
    results.sfe=meanny-meannx*results.beta-kron(et,intercept);
    results.tfe=meanty-meantx*results.beta-kron(en,intercept);
    results.yhat=x*results.beta+kron(en,results.sfe)+kron(et,results.tfe)+kron(ent,intercept);
    tnvar=nvar+N+T-1;
else
    results.yhat=x*results.beta;
    tnvar=nvar;
end    
results.resid=y-results.yhat;
yme=y-mean(y);

B = (speye(N) - rho*sparse(W));
Be=zeros(nobs,1);

for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    Be([t1:t2],1)= B*results.resid([t1:t2],1);
end
epe = Be'*Be;
rsqr1=yme'*yme;
results.rsqr=1.0-epe/rsqr1; %rsquared
rsqr3 = epe/(nobs-tnvar);
rsqr2 = rsqr1/(nobs-1.0);
results.rbar = 1 - (rsqr3/rsqr2); % rbar-squared
results.sige = (1/nobs)*epe;
results.tnvar=tnvar;

sigu = epe;
sige = results.sige;
parm = [results.beta
        results.rho
        results.sige];

results.lik = f2_sempanel(parm,ywith,xwith,W,detval,T); %Elhorst
    
if N <= 500, % t-stats using information matrix (Anselin, 1982, 1988), 
%              adjusted by Elhorst for spatial panels 

t0 = clock;
B = (speye(N) - rho*sparse(W));
BI = inv(B); WB = W*BI;
pterm = trace(WB'*WB);
xpx = zeros(nvar+2,nvar+2);
xpx(1:nvar,1:nvar) = (1/sige)*xs'*xs;
% rho, rho
xpx(nvar+1,nvar+1) = T*(trace(WB'*WB) + pterm);
% sige, sige
xpx(nvar+2,nvar+2) = nobs/(2*sige*sige);
% rho, sige
xpx(nvar+1,nvar+2) = -(T/sige)*(rho*trace(WB'*WB) - trace(BI'*WB));
xpx(nvar+2,nvar+1) = xpx(nvar+1,nvar+2);
tmp = diag(inv(xpx));
bvec = [results.beta
        results.rho];
results.tstat = bvec./(sqrt(tmp(1:nvar+1,1)));
time3 = etime(clock,t0);

elseif N > 500

t0 = clock;
hessn = hessian('f2_sempanel',parm,ywith,xwith,W,detval,T); %Elhorst

if hessn(nvar+2,nvar+2) == 0
 hessn(nvar+2,nvar+2) = 1/sige;  % this is a hack for very large models that 
end;                             % should not affect inference in these cases


xpxi = inv(-hessn); 
xpxi = diag(xpxi(1:nvar+1,1:nvar+1));
zip = find(xpxi <= 0);
 if length(zip) > 0
 xpxi(zip,1) = 1;
 fprintf(1,'sem: negative or zero variance from numerical hessian \n');
 fprintf(1,'sem: replacing t-stat with 0 \n');
 end;
 tmp = [results.beta
        results.rho];
 results.tstat = tmp./sqrt(xpxi);
 if length(zip) ~= 0
 results.tstat(zip,1) = 0;
 end;
time3 = etime(clock,t0);

end; % end of t-stat calculations

results.lndet = detval;

results.time = etime(clock,timet);
results.time1 = time1;
results.time2 = time2;
results.time3 = time3;
results.time4 = time4;

function [rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,iter,options] = sem_parse(info)
% PURPOSE: parses input arguments for far, far_g models
% ---------------------------------------------------
%  USAGE: [rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,iter] = far_parse(info)
% where info contains the structure variable with inputs 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------

% set defaults
options = optimset('fminbnd');
options.MaxIter = 500;

eflag = 1;     % default to not computing eigenvalues
ldetflag = 1;  % default to 1999 Pace and Barry MC determinant approx
order = 50;    % there are parameters used by the MC det approx
iter = 30;     % defaults based on Pace and Barry recommendation
rmin = -0.99;     % use -1,1 rho interval as default
rmax = 0.99;
detval = 0;    % just a flag
convg = 0.0001;
maxit = 500;

fields = fieldnames(info);
nf = length(fields);
if nf > 0
    
 for i=1:nf
    if strcmp(fields{i},'rmin')
        rmin = info.rmin;  eflag = 1;
    elseif strcmp(fields{i},'rmax')
        rmax = info.rmax;  eflag = 1;
    elseif strcmp(fields{i},'eigs')
        eflag = info.eigs; % flag for compute the eigenvalues
    elseif strcmp(fields{i},'convg')
       options.TolFun = info.convg;
    elseif strcmp(fields{i},'maxit')
        options.MaxIter = info.maxit;  
    elseif strcmp(fields{i},'lndet')
    detval = info.lndet;
    ldetflag = -1;
    eflag = 1;
    rmin = detval(1,1);
    nr = length(detval);
    rmax = detval(nr,1);
    elseif strcmp(fields{i},'lflag')
        tst = info.lflag;
        if tst == 0,
        ldetflag = 0; 
        elseif tst == 1,
        ldetflag = 1; 
        elseif tst == 2,
        ldetflag = 2; 
        else
        error('sar: unrecognizable lflag value on input');
        end;
    elseif strcmp(fields{i},'order')
        order = info.order;  
    elseif strcmp(fields{i},'iter')
        iter = info.iter; 
    end;
 end;
 
else, % the user has input a blank info structure
      % so we use the defaults
end; 

function [rmin,rmax,time2] = sem_eigs(eflag,W,rmin,rmax,n);
% PURPOSE: compute the eigenvalues for the weight matrix
% ---------------------------------------------------
%  USAGE: [rmin,rmax,time2] = far_eigs(eflag,W,rmin,rmax,W)
% where eflag is an input flag, W is the weight matrix
%       rmin,rmax may be used as default outputs
% and the outputs are either user-inputs or default values
% ---------------------------------------------------


if eflag == 0
t0 = clock;
opt.tol = 1e-3; opt.disp = 0;
lambda = eigs(sparse(W),speye(n),1,'SR',opt);  
rmin = 1/lambda;   
rmax = 1;
time2 = etime(clock,t0);
else
time2 = 0;
end;


function [detval,time1] = sem_lndet(ldetflag,W,rmin,rmax,detval,order,iter);
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
            error('sem: wrong lndet input argument');
        end;
        [n1,n2] = size(detval);
        if n2 ~= 2
            error('sem: wrong sized lndet input argument');
        elseif n1 == 1
            error('sem: wrong sized lndet input argument');
        end;          
end;

