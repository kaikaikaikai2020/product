function [data, Ht] = diagonal_bekk_simulate(t,k,parameters,p,q);
% PURPOSE:
%     To simulate a diagonal BEKK multivariate GARCH model. 
% 
% USAGE:
%     [data, Ht] = diagonal_bekk_simulate(t,parameters,p,q)
% 
% INPUTS:
%     t             - Length of data serie to prouce
%     k             - Dimension of series to produce
%     p             - The lag length of the innovation process
%     q             - The lag length of the AR process
%     options       - (optional) Options for the optimization(fminunc)
% 
% OUTPUTS:
%     data          - A t by k matrix of zero mean residuals
%     Ht            - A k x k x t 3 dimension matrix of conditional covariances
% 
% 
% COMMENTS:
% 
% Author: Kevin Sheppard
% kksheppard@ucsd.edu
% Revision: 2    Date: 12/31/2001
if length(parameters)~=(k*(k+1)/2)+p*k+q*k
    error('Wrong number of parameters')
end

t=t+500;

%Reshape the parameters
C=zeros(k,k);
A=zeros(k,k,p);
B=zeros(k,k,q);

C=ivech(parameters(1:(k*(k+1)/2))');
C=tril(C);
holder=(k*(k+1)/2);

LHS=eye(k);
for i=1:p
A(:,:,i)=diag(parameters(holder+1:holder+k));    
LHS=LHS-A(:,:,i)'*A(:,:,i);
holder=holder+k;
end

for i=1:p
B(:,:,i)=diag(parameters(holder+1:holder+k));    
LHS=LHS-B(:,:,i)'*B(:,:,i);
holder=holder+k;
end

Ht=zeros(k,k,t);
%Calculate the unconditional covariance
U=LHS^(-1)*(C*C')
data=randn(t,k)*U^(0.5);
m=max(p,q);
Ht(:,:,1:m)=U;

const=C*C';
for i=m+1:t+m;
    Ht(:,:,i)=const;
    for j=1:p
         Ht(:,:,i)=Ht(:,:,i)+A(:,:,j)*(data(i-j,:))'*(data(i-j,:))*A(:,:,j)';
    end
    for j=1:q
         Ht(:,:,i)=Ht(:,:,i)+B(:,:,j)*Ht(:,:,i-j)*B(:,:,j)';
    end
    data(i,:)=randn(1,k)*Ht(:,:,i)^(0.5);
end
data=data(m+500:t+m,:);
Ht=Ht(:,:,m+500:t+m);
