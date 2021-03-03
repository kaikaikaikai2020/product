%clear;
% Daily data on EWA-EWC
addpath(genpath(fullfile(pwd,'jplv7')))

%{
load('inputData_ETF', 'tday', 'syms', 'cl');
idxA=find(strcmp('EWA', syms));
idxC=find(strcmp('EWC', syms));

x=cl(:, idxA);
y=cl(:, idxC);
%}
%sql_str = 'select tradeDate,closePrice from polygon.forex_day where ticker = "%s" order by tradeDate';
sql_str = 'select tradeDate,closePrice from aksharedata.currency_hist where ticker = "%s" and tradeDate>="2011-01-01" order by tradeDate';

%usdtwd usdcnh
symbol_pool = {'usdtwd','usdkrw';'usdtwd','usdcnh';'usdtry','usdzar'};
delta_pool = [1e-4,1e-7];
re = cell(2,1);
for pool_id = 1:2
    sym1 = symbol_pool{pool_id,1};
    sym2 = symbol_pool{pool_id,2};

    x=fetchmysql(sprintf(sql_str,sym1),2);
    y=fetchmysql(sprintf(sql_str,sym2),2);
    re{pool_id} = {x,y};
end

save tempdata20200910b re