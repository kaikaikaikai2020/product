%{
当9988 HK在香港开始交易的时候（大概是2019的某天），
1） 每天在尾盘卖出。执行价= close
2)    当天晚上美国开盘的时候，以开盘5分钟的价格买入BABA US
3）取得当天的USDHKD FX的价格
4）计算pnl。1 股BABA = 8股9988 HK。 PNL = -8 x (9988的执行价格）/ (USDHKD FX) + BABA 的执行价格。

我们先开始算阿里巴巴，看看结果如何。

1 导入三组数据
2 核查数据时间点是否正确
3 公式 PNL = -8 x (9988的执行价格）/ (USDHKD FX) + BABA 计算

%}


sql_str1 = 'select tradeDate,closePrice from yuqerdata.MktHKEqudGetS54 where ticker = "09988" order by tradeDate';
x1 = fetchmysql(sql_str1,2);

%sql_str2 = 'select tradeDate,openPrice from polygon.usastock_day where ticker = "BABA" and tradeDate>="2019-11-26" order by tradeDate';
sql_str2 = 'select date(tradeDate) as t,avg(openPrice) from polygon_stock_minute.BABA where  tradeDate>="2019-11-01" and hour(tradeDate)=21 and minute(tradeDate)>=30 and minute(tradeDate)<=34 group by t;';
x2 = fetchmysql(sql_str2,2);

sql_str3 = 'select tradeDate,openPrice from polygon.forex_day where ticker = "USDHKD" and tradeDate>="2019-11-26" order by tradeDate';
x3 = fetchmysql(sql_str3,2);


[inds,commValue] = suscc_intersect({x1(:,1),x2(:,1),x3(:,1)});

x1 = x1(inds(:,1),:);
x2 = x2(inds(:,2),:);
x3 = x3(inds(:,3),:);


tref = x1(:,1);
x1 = cell2mat(x1(:,2));
x2 = cell2mat(x2(:,2));
x3 = cell2mat(x3(:,2));

pnl = -8*x1./x3+x2;
pct = -pnl./x1/8;
yc = cumprod(pct+1);
figure_S53(yc,tref,[]);
