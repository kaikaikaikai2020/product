%{
查看主要指数开首盘变化时间

%}

clear

write_sel = false;
key_str = '主要指数测试';

dN = 'S40_america';
symbols = {'mdax','aex','cac40','dax','estx50','xjo','nk200'};

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
    'hour(tradingdate),minute(tradingdate),openprice,closeprice from %s.%s ',...
    'where tradingdate>=''2010-01-01'' order by tradingdate'];
    title_str = symbols{i_sym};
    sub_sql_str = sprintf(sql_str,dN,title_str);
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
    for i = 1:T
        sub_ind = eq(day_tick,day_tick_u(i));
        sub_tt = temp_t(sub_ind,:);
        x1(i) = min(sub_tt);
        y1(i) = max(sub_tt);
    end
    figure;
    plot([x1,y1])
    title(title_str)
    
end

save('foreign_main_index.mat','re','symbols')