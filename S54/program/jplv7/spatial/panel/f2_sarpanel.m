function llike = f2_sarpanel(parm,y,x,W,detval,T)
% PURPOSE: evaluates log-likelihood -- given ML estimates
%  spatial panel autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:llike = f2_sar(parm,y,X,W,ldet,T)
%  where: parm = vector of maximum likelihood parameters
%                parm(1:k-2,1) = b, parm(k-1,1) = rho, parm(k,1) = sige
%         y    = dependent variable vector (n x 1)
%         X    = explanatory variables matrix (n x k)
%         W    = spatial weight matrix
%         ldet = matrix with [rho log determinant] values
%                computed in sar.m using one of Kelley Pace's routines  
%         T    = number of time points
% ---------------------------------------------------
%  RETURNS: a  scalar equal to minus the log-likelihood
%           function value at the ML parameters
% ---------------------------------------------------

% written by: James P. LeSage 1/2000
% University of Toledo
% Department of Economics
% Toledo, OH 43606
% jlesage@spatial.econometrics.com

% partly rewritten by J.P. Elhorst 9/2004 to account for spatial panels
% "Specification and Estimation of Spatial Panel Data Models",
% International Regional Science Review, Vol. 26, pp. 244-268.

n = length(y); 
k = length(parm);
b = parm(1:k-2,1);
rho = parm(k-1,1);
sige = parm(k,1);
n1=n/T;

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

e=y-x*b;
for t=1:T
    t1=1+(t-1)*n1;t2=t*n1;
    e([t1:t2],1)= e([t1:t2],1)-rho*sparse(W)*y([t1:t2],1);
end

epe = e'*e;
tmp2 = 1/(2*sige);
llike = -(n/2)*log(2*pi*sige) + T*detm - tmp2*epe;