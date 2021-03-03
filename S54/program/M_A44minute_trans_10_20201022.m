%{
我们有美国和中国的1分钟数据，请您帮忙测试一下
用open价格算出来仓位但执行的价格是开盘之后5分钟的平均价。
平仓价格是当天close价格。

如果不是收盘价做为平仓价，而是用5.20-5.30的平均价做为平仓价
%}
%%{
clear

tns = {'s54minind_nifty'};
t0_s = datenum(0,0,0,15,0,0);
sub_t = datenum(0,0,0,0,10,0);
sub_t0 = datenum(0,0,0,0,1,0);

tt = t0_s;
while tt(end)<=datenum(0,0,0,18,1,0)
    tt = cat(1,tt,tt(end)+sub_t);
end
%tt(tt>datenum(0,0,0,18,1,0))= [];
tt = cat(1,tt,datenum(0,0,0,19,20,0));
sql_str = ['select date(tradeDate)  as t1,ticker,avg(closePrice) from data_pro.%s where (hour(tradeDate)*100 + ',...
        'minute(tradeDate) >=%d) and (hour(tradeDate)*100 +minute(tradeDate)<=%d) group by t1,ticker'];

sql_str1 = 'select date(tradeDate)  as t1,ticker,closePrice from data_pro.%s where hour(tradeDate) = 18 and minute(tradeDate) = 00';
    
T = length(tns);
datatrans = [];
index_pool =  cell(T,1);
i=1;
%for i = 1:T
    for j = 1:length(tt)
        if j < length(tt)
            sub_sub_t1 = tt(j)+sub_t0;
            sub_sub_t2 = tt(j+1);

            v1 = hour(sub_sub_t1)*100+minute(sub_sub_t1);
            v2 = hour(sub_sub_t2)*100+minute(sub_sub_t2);

            X=fetchmysql(sprintf(sql_str,tns{i},v1,v2),2);
        else
            v1 = 1901;
            X=fetchmysql(sprintf(sql_str1,tns{i}),2);
        end
        Y=cell2table(X,'VariableNames',{'date','code','close_adj'});
        Y=unstack(Y,'close_adj','code');
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
        datatrans.(sprintf('%s%d',sub_key,v1)) = temp;
        index_pool{j} = sprintf('%s%d',sub_key,v1);
    end
%end

datatrans.index=index_pool;
save('transdata20201022.mat','datatrans');
