%获取主力合约
clear
ticker = 'NQ';
sql_str = ['select UTCDate,closeprice as close from us_stock.futureminute_m ',...
    'where Ticker like "%s%%" order by LocalDate desc'];


x = fetchmysql(sprintf(sql_str,ticker),2);
tref = unique(x(:,1));
y = zeros(size(tref));
for i = 1:length(tref)
    ind = strcmp(x(:,1),tref(i));
    tmp = cell2mat(x(ind,2));
    y(i) = tmp(1);
end

x=[tref,num2cell(y)];
x = cell2table(x,'VariableNames',{'date','close'});
writetable(x,sprintf('%s.csv',ticker));
