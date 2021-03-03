function h=garchcore2(data,parameters,stdEstimate,stdEstimate2,p,q,m,T,breakpt);
% PURPOSE:
%     Forward recursion to construct h's
% 
% USAGE:
%     h=garchcore(data,parameters,stdEstimate,p,q,m,T);
% 
% INPUTS:
%     See garchlikelihood
% 
% OUTPUTS:
%     See garchlikelihood
% 
% COMMENTS:
%     Helper function part of UCSD_GARCH toolbox. Used if you do not use the MEX file.
%     You should use the MEX file.
% 
% 
% Author: Kevin Sheppard
% kksheppard@ucsd.edu
% Revision: 2    Date: 12/31/2001
h=zeros(size(data));
h(1:m)=data(1)^2;
params1=parameters([1 3:length(parameters)]);
params2=parameters([2 3:length(parameters)]);
for t = (m + 1):breakpt+m
   h(t) = params1' * [1 ; data(t-(1:p)).^2;  h(t-(1:q)) ];
end
htemp=h(t);
temp=t;
h(t)=stdEstimate2.^2*(1-sum(parameters(3:length(parameters))));
for t = (breakpt+ m + 1):T
   h(t) = params2' * [1 ; data(t-(1:p)).^2;  h(t-(1:q)) ];
end
h(temp)=htemp;
