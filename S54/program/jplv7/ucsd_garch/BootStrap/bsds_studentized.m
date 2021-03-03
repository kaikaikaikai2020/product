function [c,u,l]=bsds_studentized(bench,models,B,w,boot)
% PURPOSE:
% Calculate Whites and Hansens p-vals for outperformance using 
% studentized residuals(possibly more powerful than standard bsds)
%
% USAGE:
% [c,u,l]=bsds(bench,models,B,w,boot)
% 
% INPUTS:
% bench  - The return series form the benchmark model
% models - The return series from each of the models used for comparrison
% B      - Number of Bootstrap replications
% w      - Desired block length
% boot   - 'STATIONARY' or 'BLOCK'.  Stationary will be used as default.
% 
% OUTPUTS:
% c - Consistent P-val(Hansen)
% u - Upper P-val(White)(Original RC P-vals)
% l - Lower P-val(Hansen)
% 
% COMMENTS:
% 
% 
% Author: Kevin Sheppard
% kksheppard@ucsd.edu
% Revision: 2    Date: 12/31/2001

[t,k]=size(models);


if size(bench,2)~=1
    error('Only 1 benchmark allowed');
end


if t~=size(models,1)
    error('Data and MOdels must have the same length');
end

if nargin==5
    if strcmp(boot,'BLOCK')
        p=1/w;
        [bsdata]=stationary_bootstrap([1:t]',1/w,B);
    else
        [bsdata]=block_bootstrap([1:t]',w,B);
    end
end

if nargin<5
    [bsdata]=block_bootstrap([1:t]',1/w,B);
end
%OK now we have teh bootstraps, what to do with them?
diffs=models-repmat(bench,1,k);
%First we need to calculate the variances of the diffs using the bootstrap
[bsdata2]=stationary_bootstrap([1:t]',p,B);
stds=zeros(1,k)
for i=1:k
    workdata=diffs(:,i);
    workdata2=mean(workdata(bsdata))';
    stds(i)=std(workdata2);
end
stat=max(mean(diffs)./stds);
% Now we need to bootstrap teh k series B times, and for each B, we need to bootstrap 
% that B times to get the variance fo that bootstrap for studentization

means=zeros(k,B);
stds2=zeros(k,B);
for i=1:B
    [bsdata3]=stationary_bootstrap([1:t]',1/w,B);
    workdata=diff(repmat(bsdata(:,i),1,k),:);
    means(i,:)=mean(workdata);
    %Now we need to bootstrap the bootstraps to get teh variance of the means
    for j=1:B
        workdata2=workdata(repmat(bsdata3(:,i),1,k),:);
        stds2(i,:)=std(workdata2);
    end
end
% Now we have B means for each series and B stds for each series coresponding to the b means

% Now we need to construct the distributions for the c,u and l





%The the consistent
%need to bootstrap the data here as well
A=zeros(k,1);
for i=1:k
    [bsdata2]=stationary_bootstrap([1:t]',1/w,B);
    temp=diffs(bsdata2,i);
    temp=mean(temp);
    A(i)=1/4*t^(0.25)*sqrt(sum((t^(0.5)*temp-t^(0.5)*mean(diffs(i,:))).^2)/B);
end

g=(mean(diffs)>-A').*mean(diffs);

%We now need to build up our matrix of means for each series and each bootstrap, k x B
perf=(means-repmat(g,B,1))./stds2;
perf=max(perf,[],2);
c=mean(perf>stat);

if nargout>1
    %Then the upper
    g=mean(diffs);
    perf=(means-repmat(g,B,1))./stds2;
    perf=max(perf,[],2);
    u=mean(perf>stat);
end


if nargout>2
    % First the lower
    g=max(0,mean(diffs));
    perf=(means-repmat(g,B,1))./stds2;
    perf=max(perf,[],2);
    l=mean(perf>stat);
end

