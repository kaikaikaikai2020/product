function llike = f2_sempanel(parm,y,x,W,detval,T)
% PURPOSE: evaluates log-likelihood -- given ML parameters
%  spatial panel error model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:llike = f2_sem(parm,y,X,W,detm,T)
%  where: parm = vector of maximum likelihood parameters
%                parm(1:k-2,1) = b, parm(k-1,1) = rho, parm(k,1) = sige
%         y    = dependent variable vector (n x 1)
%         X    = explanatory variables matrix (n x k)
%         W    = spatial weight matrix
%         ldet = matrix with [rho log determinant] values
%                computed in sem_panel.m using one of Kelley Pace's routines
%         T    = number of time points
% ---------------------------------------------------                                           
%  RETURNS: a  scalar equal to minus the log-likelihood
%           function value at the ML parameters
% ---------------------------------------------------

% written by: James P. LeSage 4/2002
% University of Toledo
% Department of Economics
% Toledo, OH 43606
% jlesage@spatial.econometrics.com

% partly rewritten by J.P. Elhorst 4/2004 to account for spatial panels
% "Specification and Estimation of Spatial Panel Data Models",
% International Regional Science Review, Vol. 26, pp. 244-268.

n = length(W);
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
B = speye(n) - rho*sparse(W);
Be=zeros(n*T,1);
for t=1:T
    t1=1+(t-1)*n;t2=t*n;
    Be([t1:t2],1)= B*(y([t1:t2],1)-x([t1:t2],:)*b);
end
epe = Be'*Be;
llike = -(n*T/2)*log(2*pi*sige) + T*detm - 1/(2*sige)*epe;