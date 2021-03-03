
%{
for i = 1:length(X)
    sub_X = X{i};
    sub_X = sub_X(:,[1,1:end]);
    sub_X(:,1) = stocks(i);
    X{i} = sub_X';
end
Y = [X{:}];

tN = 'S54.A41data';
var_info = {'ticker','tradeDate','openPrice','highPrice','lowPrice','closePrice'};
datainsert_adair(tN,var_info,Y);
%}
%insert into table1(字段) select table2：从一个表复制数据，然后把数据插入到一个已存在的表中

sprintf('开始升级A71数据')
load('inputDataOHLCDaily_stocks_20120424', 'stocks');
a = sprintf('"%s"',strjoin(stocks,'","'));
sql_str = ['select ticker,tradeDate,openPrice,highPrice,lowPrice,closePrice from polygon.usastock_day ',...
'where tradeDate> "2020-07-25" and ticker in (%s)'];
x = fetchmysql(sprintf(sql_str,a),2);
%dN = 'S54';
tN = 'S54.A41data';
var_info = {'ticker','tradeDate','openPrice','highPrice','lowPrice','closePrice'};
if ~isempty(x)
    datainsert_adair(tN,var_info,x)
end