%{
国外指数
如何划分？
有时间划分的按照早晨-下午
中间连续的，就从中间找一个节点划分

期货-手续费为十万分之五，外汇为十万分之1

{'LIGHTCMDUSD': 'WTI原油',
 'XAGUSD': 'XAGUSD',
 'XAUUSD': 'XAUUSD',
 'GBRIDXGBP': '富时100伦敦指数',
 'bulunteyuanyou': '布伦特原油',
 'dezhi30': '德指30',
 'ITAIDXEUR': '意大利富时40',
 'JPNIDXJPY': '日经225指数',
 'EUSIDXEUR': '欧洲50指数',
 'FRAIDXEUR': '法国40指数',
 'AUSIDXAUD': '澳指',
 'CHEIDXCHF': '瑞士20指数',
 'nasidake100': '纳斯达克100',
 'DOLLARIDXUSD': '美元指数',
 'NLDIDXEUR': '荷兰25指数',
 'ESPIDXEUR': '西班牙35指数'}

%}

clear

write_sel = false;
key_str = '国外商品期货测试';

dN = 'foreign_index_min_V2';
symbols = fetchmysql(sprintf('show tables from %s',dN),2);
T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);

if write_sel
    obj_wd = wordcom(fullfile(pwd,sprintf('%s.doc',key_str)));
end

tradeDay_av = zeros(T_symbols,1);
for i_sym = 1% 1:T_symbols

    P = [];
    P.feeOpen=5/100000;
    P.feeClose=5/100000;
    P.matchRecord=1;%匹配数据源：沪深300
    P.tradeRecord=1;%交易数据源：股指期货主力合约
    P.tradeMin=135;%使用早盘135分钟K线数据进行分形匹配
    P.dayMin=225;%每个交易日共225根1分钟K线
    P.M=20;%找M个最为相似的交易日
    P.muchPara=0.5;%多数上涨或下跌比例
    P.deanMethod=3;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离
    P.stopMethod=3;%1收盘价止损/2触发价止损，不考虑能否交易/3跳开则开盘价止损，否则触发价止损
    P.testStart='2010-4-16';
    P.trade_mode = 2;%1只多仓 2 多仓和空仓

    title_str = symbols{i_sym};  
    sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
        'hour(tradingdate),minute(tradingdate),open,close from %s.%s ',...
        'order by tradingdate limit 10'];
    sub_sql_str = sprintf(sql_str,dN,title_str);
    x = fetchmysql(sub_sql_str);
    if isempty(x)
        sprintf('%s %s 时间不足4年，跳出',key_str,title_str)
        continue
    end
    temp = x(:,4)*100+x(:,5);
    id_sel = temp>=900&temp<=1500 ;
    %id_sel = (temp>=800&temp<=1000) |(temp>=1100&temp<=1430);
    x = x(id_sel,:);
    if isempty(x)
        continue
    end
    t = datenum([x(:,1:5),zeros(size(x(:,1)))]);
    %统计中间有停牌的情况，并剔除
    day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
    day_tick_u = unique(day_tick);
    ind_miss = false(size(day_tick));
    ind_miss_u = false(size(day_tick_u));
    T = length(day_tick_u);
    temp_tradedays = zeros(T,1);
    for i = 1:T
        sub_ind = eq(day_tick,day_tick_u(i));
        if sum(sub_ind)<225
            ind_miss(sub_ind) = true;
            ind_miss_u(i) = true;
        end
        temp_tradedays(i) = sum(sub_ind);
    end
    temp_tradedays(ind_miss_u) = [];
    if ~isempty(temp_tradedays)
        tradeDay_av(i_sym) = mean(temp_tradedays);
    end
    
    day_tick_u(ind_miss_u) = [];
    x(ind_miss,:) = [];
    t(ind_miss,:) = [];

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
    [tradeYield,result1,tradeDetail,yearDetail,h] = SMTTradingModelTool(closeprice,openprice,t,closeprice,openprice,t,P);
    ah = gca;
    title(ah,title_str2);
    if write_sel
        obj_wd.pasteFigure(h,title_str2);
    end
    y_c = cumprod(tradeYield(:,2)+1);
    %统计参数
    [v0,v_str0] = curve_static(y_c,[],false);
    [v,v_str] = ad_trans_sta_info(v0,v_str0); 
    result2 = [v_str;v]';
    result = [{'',title_str2};[result1;result2]];
    sta_re{i_sym} = result;
    sprintf('%s %d-%d',key_str,i_sym,T_symbols)
    
end
if write_sel
    obj_wd.CloseWord()
end
y = [sta_re{:}];
y = y(:,[1,2:2:end]);
y = y';
