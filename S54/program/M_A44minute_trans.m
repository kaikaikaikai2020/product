%{
我们有美国和中国的1分钟数据，请您帮忙测试一下
用open价格算出来仓位但执行的价格是开盘之后5分钟的平均价。
平仓价格是当天close价格。
%}
%%{
clear

tns = {'s54minind_as51';'s54minind_bse100';'s54minind_hsi';'s54minind_krx';'s54minind_msci';'s54minind_nky';'s54minind_sx5e';'s54minind_hscei';'s54minind_nifty'};
t0_s=[8,0;11,37;9,19;8,0;9,0;8,0;15,0;9,19;11,37];
sql_str = ['select date(tradeDate)  as t1,ticker,avg(openPrice) from data_pro.%s where hour(tradeDate)=%d ',...
        'and minute(tradeDate)>=%d and minute(tradeDate)<=%d group by t1,ticker'];

T = length(tns);
datatrans = [];
index_pool =  cell(T,1);
for i = 1:T
    sub_t0_s=t0_s;
    sub_h = sub_t0_s(i,1);
    sub_t=sub_t0_s(i,2);
    X=fetchmysql(sprintf(sql_str,tns{i},sub_h,sub_t,sub_t+4),2);
    
    Y=cell2table(X,'VariableNames',{'date','code','open_adj'});
    Y=unstack(Y,'open_adj','code');
    tref=Y(:,1);
    stocks = Y.Properties.VariableNames(2:end);
    X=table2array(Y(:,2:end));
    tref = table2cell(tref);
    sprintf('%d-%d',i,T)
    
    sub_key =strsplit( tns{i},'_');
    sub_key = upper(sub_key{2});
    temp=[];
    temp.tref = tref;
    temp.stocks = stocks;
    temp.X = X;
    datatrans.(sub_key) = temp;
    index_pool{i} = sub_key;
end

datatrans.index=index_pool;
save('transdata20201011.mat','datatrans');
