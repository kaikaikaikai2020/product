%{
因为2005-2009的数据其实不是很重要，如果是后面该了时间，我们能不能动态的检测开始时间？
因为这个欧洲时间可能又夏令时冬令时的问题。 
对于cac/dax/estx这个可不可以自动检测开盘时间，
取前4个小时比较距离，后四个小时执行试一下

%}

clear

write_sel = true;
key_str = '主要指数测试-固定时间点';

dN = 'S40_america';
symbols = {'cac40','dax','estx50'};
T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);

if write_sel
    obj_wd = wordcom(fullfile(pwd,sprintf('%s.doc',key_str)));
end

tradeDay_av = zeros(T_symbols,1);
re = cell(T_symbols,1);

min_T_num = 60*4;
for i_sym =  1:T_symbols
    P = [];
    P.feeOpen=5/100000;
    P.feeClose=5/100000;
    P.matchRecord=1;%匹配数据源：沪深300
    P.tradeRecord=1;%交易数据源：股指期货主力合约
    P.tradeMin=min_T_num;%使用早盘135分钟K线数据进行分形匹配
    P.dayMin=min_T_num*2;%每个交易日共225根1分钟K线
    P.M=20;%找M个最为相似的交易日
    P.muchPara=0.5;%多数上涨或下跌比例
    P.deanMethod=3;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离
    P.stopMethod=3;%1收盘价止损/2触发价止损，不考虑能否交易/3跳开则开盘价止损，否则触发价止损
    P.testStart='2010-4-16';
    P.trade_mode = 2;%1只多仓 2 多仓和空仓

    title_str = symbols{i_sym};  
    if ~strcmpi(title_str,'nk200')
        sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
        'hour(tradingdate),minute(tradingdate),openprice,closeprice from %s.%s ',...
        ' order by tradingdate'];
    else
        sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
        'hour(tradingdate),minute(tradingdate),openprice,closeprice from %s.%s ',...
        'where tradingdate>=''2010-01-01'' order by tradingdate'];
    end
    sub_sql_str = sprintf(sql_str,dN,title_str);
    x = fetchmysql(sub_sql_str);
    if isempty(x)
        sprintf('%s %s 时间不足4年，跳出',key_str,title_str)
        continue
    end
    temp_t = x(:,4)*100+x(:,5);
    if isempty(x)
        continue
    end
    %t = datenum([x(:,1:5),zeros(size(x(:,1)))]);
    %统计中间有停牌的情况，并剔除
    day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
    day_tick_u = unique(day_tick);
    ind_miss = false(size(day_tick));
    ind_miss_u = false(size(day_tick_u));
    T = length(day_tick_u);
    temp_tradedays = zeros(T,1);
    x1 = cell(T,1);
    for i = 1:T
        sub_ind = eq(day_tick,day_tick_u(i));
        sub_x = x(sub_ind,:);
        sub_tt = temp_t(sub_ind,:);
        if sum(sub_ind)<min_T_num*2
            ind_miss(sub_ind) = true;
            ind_miss_u(i) = true;
        else
           %filling 使用前面的值
            sub_sub_x1 = sub_x([1:min_T_num,end-min_T_num+1:end],:);
            x1{i} = sub_sub_x1';
        end
        temp_tradedays(i) = sum(sub_ind);
    end
    temp_tradedays(ind_miss_u) = [];
    if ~isempty(temp_tradedays)
        tradeDay_av(i_sym) = mean(temp_tradedays);
    end
    x1(ind_miss_u) = [];
    
    day_tick_u(ind_miss_u) = [];
    x = [x1{:}]';
    if isempty(x)
        continue
    end
    t = datenum([x(:,1:5),zeros(size(x(:,1)))]);

    %初始时间设定
    min_day_num = 210*4;
    if length(day_tick_u)<min_day_num
        sprintf('%s %s 时间不足4年，跳出',key_str,title_str)
        continue
    else
        temp = num2str(day_tick_u(min_day_num/2+1));
        P.testStart=sprintf('%s-%s-%s',temp(1:4),temp(5:6),temp(7:8));
    end

    openprice = x(:,end-1);
    closeprice = x(:,end);
    %SMTTradingModelTool(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate)
    [tradeYield,result1,tradeDetail,yearDetail,h,assure_ratio] = SMTTradingModelTool(closeprice,openprice,t,closeprice,openprice,t,P);
    re{i_sym} = [tradeYield,assure_ratio];
    ah = gca;
    title(ah,title_str);
    if write_sel
        obj_wd.pasteFigure(h,title_str);
    end
    y_c = cumprod(tradeYield(:,2)+1);
    %统计参数
    [v0,v_str0] = curve_static(y_c,[],false);
    [v,v_str] = ad_trans_sta_info(v0,v_str0); 
    result2 = [v_str;v]';
    result = [{'',title_str};[result1;result2]];
    sta_re{i_sym} = result;
    sprintf('%s %d-%d',key_str,i_sym,T_symbols)
    
end
if write_sel
    obj_wd.CloseWord()
end
y = [sta_re{:}];
y = y(:,[1,2:2:end]);
y = y';
xlswrite(sprintf('%s.xlsx',key_str),y);

save('foreign_main_index.mat','re','symbols')