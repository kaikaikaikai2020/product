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

1���۹�β������8��
2�����ɿ�������1��
3��������������1��
4���۹ɿ�������8��

%}
par_info = {"01179","HTHT";"09618","JD";"06160","BGNE";"09991","BZUN";"09688","ZLAB"};
forex_pool = {'USDHKD','USDHKD','USDHKD','USDHKD','USDHKD'};
r = [1,2,13,3,1];

% data = load('NTES_20201018.mat');
% data = data.sub_data;
% for i = 1:length(data.stocks)
%     temp = data.stocks{i};
%     temp = split(temp,' ');
%     data.stocks{i} = temp{1};
% end

for isel =  1:size(par_info,1)
    ticker1 = par_info{isel,1};
    ticker2 = par_info{isel,2};
    sub_r = r(isel);
    sub_forex=forex_pool{isel};
    sql_str1 = 'select tradeDate,openPrice,closePrice from yuqerdata.MktHKEqudGetS54 where ticker = "%s" order by tradeDate';
    x1 = fetchmysql(sprintf(sql_str1,ticker1),2);
%     if eq(isel,1) || eq(isel,2)
%         %sql_str2 = 'select tradeDate,openPrice,closePrice from yuqerdata.MktUsequdGetS54 where ticker = "%s" and tradeDate>="2019-11-26" order by tradeDate';
%         %polygon.usastock_day
    sql_str2 = 'select tradeDate,openPrice,closePrice from yuqerdata.MktUsequdGetS54 where ticker = "%s" and tradeDate>="2019-11-26" order by tradeDate';
    %sql_str2 = 'select date(tradeDate) as t,avg(openPrice) from polygon_stock_minute.BABA where  tradeDate>="2019-11-01" and hour(tradeDate)=21 and minute(tradeDate)>=30 and minute(tradeDate)<=34 group by t;';
    x2 = fetchmysql(sprintf(sql_str2,ticker2),2);
%     %elseif eq(isel,2)
%     %    ind = strcmp(data.stocks,ticker2);
%     %    x2=[data.tref,num2cell([data.PX_OPEN(:,ind),data.PX_LAST(:,ind)])];
%     elseif eq(isel,3)
%         temp = load(ticker2);
%         x2 = [temp.X.tradeDate,num2cell([temp.X.openPrice,temp.X.closePrice])];
%     end
    sql_str3 = 'select tradeDate,openPrice from polygon.forex_day where ticker = "%s" and tradeDate>="2019-11-26" order by tradeDate';
    x3 = fetchmysql(sprintf(sql_str3,sub_forex),2);
    [inds,commValue] = suscc_intersect({x1(:,1),x2(:,1),x3(:,1)});

    x1 = x1(inds(:,1),:);
    x2 = x2(inds(:,2),:);
    x3 = x3(inds(:,3),:);


    tref = x1(:,1);
    x1 = cell2mat(x1(:,2:end));
    x2 = cell2mat(x2(:,2:end));
    x3 = cell2mat(x3(:,2:end));

    pnl = -sub_r*x1(:,2)./x3+x2(:,1);
    pct1 = -pnl./x1(:,2)/sub_r;

    pnl2 = x2(1:end-1,2)-sub_r*x1(2:end,1)./x3(2:end);
    pct2 = -pnl2./x2(1:end-1,2);
    pct2 = [0;pct2];
    pct = (pct1+pct2)/2;
    y_re =cumprod(1+[pct1,pct2,pct]);
    h = figure_S53(y_re,tref,[]);
    legend({'�۹��������ɿ���','�������̸۹ɿ���','�ۺ�'},'location','best')


    re = [tref,num2cell([x1,x2,x3])];
    re = cell2table(re,'VariableNames',{'tradeDate','HKopenPrice','HKclosePrice','USopenPrice','USclosePrice','USDHKD'});
    writetable(re,sprintf('position%s%s.csv',ticker1,ticker2));
    
    title(sprintf('%s-%s',ticker1,ticker2))
end


