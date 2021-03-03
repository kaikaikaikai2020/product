function [LLF,likelihoods,Ht]=diagonal_bekk_likelihood(parameters,errors,p,q,k,k2,t);
% PURPOSE:
%      To Estimate a full BEKK multivariate GARCH likelihood.
% 
% USAGE:
%      [LLF,likelihoods,Ht]=full_bekk_mvgarch_likelihood(parameters,errors,p,q,k,k2,t);
% 
% INPUTS:
%      parameters - a k*(k+1)/2 + k^2 *p +k^2*q vector of model parameters of the form
%                   [ivech(C);reshape(A(1),k*k,1);...;reshape(A(p),k*k,1);reshape(B,k*k,1); ...
%                     reshape(B(q),k*k,1)]
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



%The first k(k+1)/2 parameters are C, the next p are A, and the next q are B
C=parameters(1:(k2));
A=parameters(k2+1:k2+k*k*p);
B=parameters(k2+k*k*p+1:k2+k*k*p+k*k*q);
tempA=zeros(k,k,p);
tempB=zeros(k,k,p);
for i=1:p
    tempA(:,:,i)=reshape(A((k*k*(i-1)+1):(k*k*i)),k,k);
end
for i=1:q
    tempB(:,:,i)=reshape(B((k*k*(i-1)+1):(k*k*i)),k,k);
end
A=tempA;
B=tempB;

C=ivech(C);
C=tril(C);
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
for i=m+1:t+m;
    Ht(:,:,i)=const;
    for j=1:p
         Ht(:,:,i)=Ht(:,:,i)+A(:,:,j)*(errors(i-j,:))'*(errors(i-j,:))*A(:,:,j)';
    end
    for j=1:q
         Ht(:,:,i)=Ht(:,:,i)+B(:,:,j)*Ht(:,:,i-j)*B(:,:,j)';
    end
    likelihoods(i)=k*log(2*pi)+(log(det(Ht(:,:,i)))+errors(i,:)*Ht(:,:,i)^(-1)*errors(i,:)');
    LLF=LLF+likelihoods(i);
end
LLF=0.5*(LLF);
likelihoods=0.5*likelihoods(m+1:t+m);
Ht=Ht(:,:,m+1:t+m);
if isnan(LLF)
    LLF=1e6;
end
