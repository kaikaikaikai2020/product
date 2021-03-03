%{
�������������й���1�������ݣ�������æ����һ��
��open�۸��������λ��ִ�еļ۸��ǿ���֮��5���ӵ�ƽ���ۡ�
ƽ�ּ۸��ǵ���close�۸�
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