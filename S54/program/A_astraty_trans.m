%{
��9988 HK����ۿ�ʼ���׵�ʱ�򣨴����2019��ĳ�죩��
1�� ÿ����β��������ִ�м�= close
2)    ���������������̵�ʱ���Կ���5���ӵļ۸�����BABA US
3��ȡ�õ����USDHKD FX�ļ۸�
4������pnl��1 ��BABA = 8��9988 HK�� PNL = -8 x (9988��ִ�м۸�/ (USDHKD FX) + BABA ��ִ�м۸�

�����ȿ�ʼ�㰢��Ͱͣ����������Ρ�

1 ������������
2 �˲�����ʱ����Ƿ���ȷ
3 ��ʽ PNL = -8 x (9988��ִ�м۸�/ (USDHKD FX) + BABA ����

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
