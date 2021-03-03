function [inds,commValue] = suscc_intersect(varargin)
x = varargin;
if eq(length(x),1)
    x = x{1};
end
commValue = x{1};
for i = 2:length(x)
    commValue = intersect(x{i},commValue);
end
inds = zeros(length(commValue),length(x));
for i = 1:length(x)
    [~,inds(:,i)] = intersect(x{i},commValue);
end