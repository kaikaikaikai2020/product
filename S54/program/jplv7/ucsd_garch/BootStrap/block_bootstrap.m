function [bsdata, indexes]=block_bootstrap(data,w,B);
% PURPOSE:
%   Implements a block bootstrap for bootstrapping stationary, dependant series
% 
% USAGE:
%     [bsdata, indexes]=block_bootstrap(data,w,B);
% 
% INPUTS:
%     data: T by 1 matrix of data to be bootstrapped(shoudl be unit root)
%     w:    Desired average windows length
%     B:    Number fo bootstraps
% 
% OUTPUTS:
%     bsdata: T x B matrix of bootstrapped data
%     indexes: T by B matrix of location sof the original data(data(indexes)=bsdata;
% 
% COMMENTS:
% 
% Author: Kevin Sheppard
% kksheppard@ucsd.edu
% Revision: 2    Date: 12/31/2001

[t,k]=size(data);
data=[data;data(1:w)];
s=ceil(t/w);
Bs=rand(s,B);
indexes=zeros(t,B);
index=1;
for i=1:t
    if mod(i,w)==1
        indexes(i,:)=ceil(Bs(index,:)*t);
        index=index+1;    
    else
        indexes(i,:)=indexes(i-1,:)+1;
    end
end

bsdata=data(indexes);


    








