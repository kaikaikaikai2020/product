clear
dN = 'S54';
tn = 'SHARADAR_SF1';
[~,~,x] = xlsread('temp.csv');
var_info = x(:,1);
var_type = cell(size(var_info));
var_type(:) = {'float'};
var_type(1:2) = {'varchar(12)'};
var_type(3:6) = {'date'};
key_var = [];
[OK1,OK2,OK3] = create_table_adair(dN,tn,var_info,var_type,key_var);

dN = 'S54';
tn = 'A71data';
var_info = {'tradeDate','closePrice','openPrice','highestPrice','lowestPrice','volume','CHG'};
var_type = cell(size(var_info));
var_type(:) = {'float'};
var_type([1,6]) = {'date','varchar(20)'};
key_var = 'tradeDate';
[OK1,OK2,OK3] = create_table_adair(dN,tn,var_info,var_type,key_var);

dN = 'S54';
tn = 'A41data';
var_info = {'ticker','tradeDate','openPrice','highPrice','lowPrice','closePrice'};
var_type = cell(size(var_info));
var_type(:) = {'float'};
var_type([1,2]) = {'varchar(12)','date'};
key_var = 'ticker,tradeDate';
[OK1,OK2,OK3] = create_table_adair(dN,tn,var_info,var_type,key_var);


