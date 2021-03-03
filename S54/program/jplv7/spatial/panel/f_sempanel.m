function lik = f_sempanel(rho,eD,W,detval,T)
% PURPOSE: evaluates concentrated log-likelihood for the 
%  spatial panel error model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:llike = f_sem(rho,eD,W,detm,T)
%  where: rho  = spatial error parameter
%         eD   = begls residuals
%         W    = spatial weight matrix
%         detm =  matrix with [rho log determinant] values
%                computed in sem_panel.m using one of 
%                Kelley Pace's routines
%         T    = number of time points
% ---------------------------------------------------
%  RETURNS: a  scalar equal to minus the log-likelihood
%           function value at the parameter rho
% ---------------------------------------------------

% written by: James P. LeSage 1/2000
% University of Toledo
% Department of Economics
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

% partly rewritten by J.P. Elhorst 4/2004 to account for spatial panels
% "Specification and Estimation of Spatial Panel Data Models",
% International Regional Science Review, Vol. 26, pp. 244-268.

n = length(W); 
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
    Be([t1:t2],1)= B*eD([t1:t2],1);
end
epe = Be'*Be;
%lik = (n*T/2)*log(2*pi*epe) - T*detm;
lik = (n*T/2)*log(epe) - T*detm;