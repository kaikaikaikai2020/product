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
close all
par_info = {"09988","BABA";"09999","NTES";"00700","TCEHY";"09618","JD"};
forex_pool = {'USDHKD','USDHKD','USDHKD','USDHKD'};
r = [8,5,1,2];

data = load('NTES_20201018.mat');
data = data.sub_data;
for i = 1:length(data.stocks)
    temp = data.stocks{i};
    temp = split(temp,' ');
    data.stocks{i} = temp{1};
end

re = cell(size(forex_pool));
temp_tref0 = load('temp_tref0.mat');
temp_tref0 = temp_tref0.tref0;
tref0 = [];
for isel =  1:size(par_info,1)
    ticker1 = par_info{isel,1};
    ticker2 = par_info{isel,2};
    sub_r = r(isel);
    sub_forex=forex_pool{isel};
    sql_str1 = 'select tradeDate,openPrice,closePrice from yuqerdata.MktHKEqudGetS54 where ticker = "%s" order by tradeDate';
    x1 = fetchmysql(sprintf(sql_str1,ticker1),2);
    if eq(isel,1) ||eq(isel,4)
        %sql_str2 = 'select tradeDate,openPrice,closePrice from yuqerdata.MktUsequdGetS54 where ticker = "%s" and tradeDate>="2019-11-26" order by tradeDate';
        %polygon.usastock_day
        sql_str2 = 'select tradeDate,openPrice,closePrice from polygon.usastock_day where ticker = "%s" and tradeDate>="2019-11-26" order by tradeDate';
        %sql_str2 = 'select date(tradeDate) as t,avg(openPrice) from polygon_stock_minute.BABA where  tradeDate>="2019-11-01" and hour(tradeDate)=21 and minute(tradeDate)>=30 and minute(tradeDate)<=34 group by t;';
        x2 = fetchmysql(sprintf(sql_str2,ticker2),2);
    elseif eq(isel,2)
       ind = strcmp(data.stocks,ticker2);
       x2=[data.tref,num2cell([data.PX_OPEN(:,ind),data.PX_LAST(:,ind)])];
    elseif eq(isel,3)
        temp = load(ticker2);
        x2 = [temp.X.tradeDate,num2cell([temp.X.openPrice,temp.X.closePrice])];
    end
    sql_str3 = 'select tradeDate,openPrice from polygon.forex_day where ticker = "%s" and tradeDate>="2019-11-26" order by tradeDate';
    x3 = fetchmysql(sprintf(sql_str3,sub_forex),2);
    [inds,commValue] = suscc_intersect({x1(:,1),x2(:,1),x3(:,1),temp_tref0});

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


    sub_re = [tref,num2cell([x1,x2,x3])];
    sub_re = cell2table(sub_re,'VariableNames',{'tradeDate','HKopenPrice','HKclosePrice','USopenPrice','USclosePrice','USDHKD'});
    writetable(sub_re,sprintf('position%s%s.csv',ticker1,ticker2));
    title(sprintf('%s-%s',ticker1,ticker2))
    
    
    sub_re1 = [tref,num2cell([pct1,pct2,pct])];
    re{isel} = sub_re1;
    if isempty(tref0)
        tref0 = sub_re1(:,1);
    else
        tref0 = intersect(tref0,sub_re1(:,1));
    end
end
re(cellfun(@isempty,re)) = [];
n = length(re);
x=zeros(length(tref0),3);
for i = 1:n
    temp = re{i};
    [~,~,ib] = intersect(tref0,temp(:,1),'stable');
    temp = cell2mat(temp(ib,2:end));
    x = x+temp;
end
x = x./n;

y_re =cumprod(1+x);
h = figure_S53(y_re,tref0,[]);
legend({'�۹��������ɿ���-��Ȩ','�������̸۹ɿ���-��Ȩ','�ۺ�-��Ȩ'},'location','best')



