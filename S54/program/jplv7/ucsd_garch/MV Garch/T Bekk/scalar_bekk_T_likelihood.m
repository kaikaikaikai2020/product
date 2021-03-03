function [LLF,likelihoods,Ht]=scalar_bekk_T_mvgarch_likelihood(parameters,errors,p,q,k,k2,t);
% PURPOSE:
%      To Estimate a scalar BEKK multivariate GARCH likelihood with T-dist errors.
% 
% USAGE:
%      [LLF,likelihoods,Ht]=scalar_bekk_T_mvgarch_likelihood(parameters,errors,p,q,k,k2,t);
% 
% INPUTS:
%      parameters - a k*(k+1)/2 + p +q +1 vector of model parameters of the form
%                   [ivech(C);(A(1));...;(A(p));(B(1)); ...(B(q));nu]
%      errors     - A zeromean t by k martix of residuals
%      p          - The lag length of the innovation process
%      q          - The lag length of the AR process
%      k          - The number of data series
%      k2         - k*k
%      t          - the length of the data series
% 
% OUTPUTS:
%      LLF           - The loglikelihood of the function at the optimum
%      Ht            - A k x k x t 3 dimension matrix of conditional covariances
%      likelihoods   - A t by 1 vector of individual likelihoods
% 
% COMMENTS:
% 
% 
% Author: Kevin Sheppard
% kksheppard@ucsd.edu
% Revision: 2    Date: 12/31/2001

 


%   Copyright (c) Kevin Sheppard.
%   $Revision: 1.0 $  $Date: 2001/04/19 $


%The first k(k+1)/2 parameters are C, the next p are A, and the next q are B
C=parameters(1:(k2));
A=parameters(k2+1:k2+p);
B=parameters(k2+p+1:k2+p+q);
nu=parameters(k2+p+q+1);
C=ivech(C);
const=C*C';

uncond=cov(errors);
% for starting up, both ee' and H have expectation uncond.  We cna leverage thsi to help the loops.
m=max(p,q);
eeprime=zeros(k,k,t+m);
Ht=zeros(k,k,t+m);
for i=1:m
    eeprime(:,:,i)=uncond;
    Ht(:,:,i)=uncond;
end
LLF=0;
errors=[repmat(diag(uncond)',m,1);errors];
likelihoods=zeros(t+m,1);

constant=(k/2)*(log(nu)-log(nu-2))+gammaln(0.5*(k+nu))-(k/2)*log(nu*pi)-gammaln(0.5*nu);
for i=m+1:t+m;
    Ht(:,:,i)=const;
    for j=1:p
         Ht(:,:,i)=Ht(:,:,i)+A(j)*(errors(i-j,:))'*(errors(i-j,:))*A(j);
    end
    for j=1:q
         Ht(:,:,i)=Ht(:,:,i)+B(j)*Ht(:,:,i-j)*B(j);
    end
    likelihoods(i)=constant-0.5*log(det(Ht(:,:,i))) - 0.5*(k+nu)*log(1+((1/(nu-2))*errors(i,:)*Ht(:,:,i)^(-1)*errors(i,:)'));
%   likelihoods(i)=k*log(2*pi)+(log(det(Ht(:,:,i)))+errors(i,:)*Ht(:,:,i)^(-1)*errors(i,:)');
    LLF=LLF+likelihoods(i);
end
LLF=-(LLF);
likelihoods=-likelihoods(m+1:t+m);
Ht=Ht(:,:,m+1:t+m);
if isnan(LLF)
    LLF=1e6;
end
if LLF < -6000
    C=parameters(1:(k2));
A=parameters(k2+1:k2+p)
B=parameters(k2+p+1:k2+p+q)
nu=parameters(k2+p+q+1)
C=ivech(C);
const=C*C'
end

