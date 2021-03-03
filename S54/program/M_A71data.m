clear
tN = 'S54.A71data';
%tn = '';
var_info = {'tradeDate','closePrice','openPrice','highestPrice','lowestPrice','volume','CHG'};
[~,~,data]= xlsread('dataA71.xlsx');
x =data(2:end,:);
x(:,1) = cellstr(datestr(datenum(x(:,1)),'yyyy-mm-dd'));
datainsert_adair(tN,var_info,x);