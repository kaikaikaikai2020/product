function llike = f_sdm(rho,y,x,W,detval)
% PURPOSE: evaluates concentrated log-likelihood for the 
%          spatial durbin model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:llike = f_sdm(rho,y,x,W,detm)
%  where: rho  = spatial autoregressive parameter
%          y   = dependent variable vector
%          x   = data matrix
%          W   = spatial weight matrix
%         detm =  matrix with [rho log determinant] values
%                computed in sdm.m using one of 
%                Kelley Pace's routines  
% ---------------------------------------------------
%  RETURNS: a  scalar equal to minus the log-likelihood
%           function value at the parameter rho
%  --------------------------------------------------
%  NOTE: this is really two functions depending
%        on nargin = 4 or nargin = 5 (see the function)
% --------------------------------------------------- 
%  SEE ALSO: sdm, f_far, f_sac, f_sem
% ---------------------------------------------------

% written by: James P. LeSage 1/2000
% University of Toledo
% Department of Economics
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

[n k] = size(x); 
rho2 = rho*rho;


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


dy=W*y;
xdx=[ x(:,2:k) W*x(:,2:k) ones(n,1)];
xdxtxdx=(xdx'*xdx);
xdxinv=inv(xdxtxdx);
xdxy=xdx'*y;
xdxdy=xdx'*dy;
bmat=xdxtxdx\[xdxy xdxdy];
bols=bmat(:,1);
bolsd=bmat(:,2);
eo=y-xdx*bols;
ed=dy-xdx*bolsd;
e2o=(eo'*eo);
edo=(ed'*eo);
e2d=(ed'*ed);
logsse=log(e2o-2*rho*edo+rho2*e2d);

llike = (n/2)*log(pi) -detm + (n/2)*logsse;
