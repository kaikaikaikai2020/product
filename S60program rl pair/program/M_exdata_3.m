clear
[~,~,x] = xlsread('pair data for RL trading.xlsx');

x([1:3,5:7],:)=[];
tref = x(2:end,1);
tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));
x = x(:,2:end);
x(:,cellfun(@isnumeric,x(1,:))) = [];
symbol_pool = x(1,:);
X = cell2mat(x(2:end,:));
T = length(symbol_pool);
for i = 1:T
    fn = symbol_pool{i};
    x = [tref,num2cell(X(:,i))];
    x = cell2table(x,'VariableNames',{'date','close'});
    fn = strsplit(fn,' ');
    fn = fn{1};
    writetable(x,sprintf('%s.csv',fn));
end
