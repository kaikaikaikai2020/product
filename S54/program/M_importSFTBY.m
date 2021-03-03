x = readtable('SFTBY历史数据.csv');
var=x.Properties.VariableNames;
x = table2cell(x);
tref = cellstr(datestr(datenum(x(:,1)),'yyyy-mm-dd'));
X=[];
X.(var{1}) = tref;
for i = 2:length(var)-1
    X.(var{i}) = cellfun(@str2double,x(:,i));
end

save SFTBY.mat X