function [bsdata, indexes]=stationary_bootstrap(data,w,B);
% PURPOSE:
%     Implements Politis' continuous bootstrap for bootstrapping unit root series
% 
% USAGE:
%     [bsdata, indexes]=stationary_bootstrap(data,w,B);
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


p=1/w;
[t,k]=size(data);
bsdata=zeros(t,B);
indexes=zeros(t,B);
indexes(1,:)=ceil(t*rand(1,B));
bsdata(1,:)=data(indexes(1,:))';

for i=2:t
    select=rand(1,B)<p;
    indexes(i,:)=ceil((t-1)*rand(1,B)).*select+(mod(~select.*(indexes(i-1,:)),t)+1);
    bsdata(i,:)=data(indexes(i,:))';
end

