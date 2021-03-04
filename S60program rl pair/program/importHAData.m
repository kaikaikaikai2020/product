function [stockList,stockData] = importHAData()
load stockList;
x = stockList(:,2);
T = length(x);
info = cell(T,1);
for i = 1:T
    temp = x{i};
    temp = split(temp,' ');
    temp = temp{1};
    info{i} = sprintf('%0.5d',str2double(temp));
end
tickerH = info;

x = stockList(:,3);
T = length(x);
info = cell(T,1);
for i = 1:T
    temp = x{i};
    temp = split(temp,' ');
    temp = temp{1};
    info{i} =temp;
end
tickerA = info;

sql_str1 = 'select tradeDate,closePrice from yuqerdata.yq_mktequdadjafget where ticker = "%s" order by tradeDate';
sql_str2 = 'select tradeDate,closePrice from data_pro.hks58 where ticker = "%s" order by tradeDate';

stockData = cell(T,1);
for i = 1:T
    x1 = fetchmysql(sprintf(sql_str2,tickerH{i}),2);
    x2 = fetchmysql(sprintf(sql_str1,tickerA{i}),2);
    
    [~,ia,ib] = intersect(x1(:,1),x2(:,1));
    x = [x1(ia,:),x2(ib,2)];
    x = [datenum(x(:,1)),cell2mat(x(:,2:end))];
    stockData{i} = x;
end