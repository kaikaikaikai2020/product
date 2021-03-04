%{
查看主要指数开首盘变化时间

%}

clear

write_sel = true;
key_str = '主要指数测试';

dN = 'S40_america';
symbols = {'HPQ';'KFT';'GE';'JNJ';'PG';'AXP';'KO';'GS';'HD';'MSFT';'XOM';...
    'TRV';'DIS';'MO';'C';'AAPL';'PFE';'MMM';'CSCO';'DOW';'AIG';'MCD';'BAC';...
    'CAT';'NKE';'V';'HON';'VZ';'AA';'CVX';'GM';'WMT';'T';'MRK';'UTX';'UNH';...
    'WBA';'BA';'DWDP';'INTC';'IBM';'JPM';'DD'};

T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);

if write_sel
    obj_wd = wordcom(fullfile(pwd,sprintf('%s.doc',key_str)));
end

tradeDay_av = zeros(T_symbols,1);
re = cell(T_symbols,1);
for i_sym =  1:T_symbols
    

    sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
    'hour(tradingdate),minute(tradingdate),openprice,closeprice from index_comp_price_temp.american_djia ',...
    'where Ticker = ''%s'' and tradingdate>=''2010-01-01'' order by tradingdate'];
    title_str = symbols{i_sym};
    sub_sql_str = sprintf(sql_str,title_str);
    x = fetchmysql(sub_sql_str);
    if isempty(x)
        sprintf('%s %s 时间不足4年，跳出',key_str,title_str)
        continue
    end
    temp_t = x(:,4)*100+x(:,5);
    %统计中间有停牌的情况，并剔除
    day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
    day_tick_u = unique(day_tick);
    ind_miss = false(size(day_tick));
    ind_miss_u = false(size(day_tick_u));
    T = length(day_tick_u);
    temp_tradedays = zeros(T,1);
    x1 = zeros(T,1);
    y1 = zeros(T,1);
    parfor i = 1:T
        sub_ind = eq(day_tick,day_tick_u(i));
        sub_tt = temp_t(sub_ind,:);
        x1(i) = min(sub_tt);
        y1(i) = max(sub_tt);
    end
    h=figure;
    plot([x1,y1])
    title(title_str)
    if write_sel
        obj_wd.pasteFigure(h,title_str);
    end
end

