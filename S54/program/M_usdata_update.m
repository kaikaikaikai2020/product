%{
我们有美国和中国的1分钟数据，请您帮忙测试一下
用open价格算出来仓位但执行的价格是开盘之后5分钟的平均价。
平仓价格是当天close价格。
%}
%%{
clear
load a41stocks.mat stocks
sql_str = ['select date(tradeDate)  as t1,ticker,avg(openPrice) from polygon_stock_minute.%s where hour(tradeDate)=13 ',...
        'and minute(tradeDate)>=30 and minute(tradeDate)<=34 group by t1'];

T = length(stocks);
X = cell(T,1);
parfor i = 1:T
    x=fetchmysql(sprintf(sql_str,stocks{i}),2);
    X{i}=x';
    sprintf('%d-%d',i,T)
end
Y=[X{:}]';
X0=X;
%Y(:,2)= cellfun(@(x) x(3:end),Y(:,2),'UniformOutput',false);

Y=cell2table(Y,'VariableNames',{'date','code','open_adj'});
%}
%X1=readtable('csidata_adj.csv');
Y=unstack(Y,'open_adj','code');
tref=Y(:,1);
stocks = Y.Properties.VariableNames(2:end);
X=table2array(Y(:,2:end));


save US_adjdata X tref stocks