%{
主要指数测试
1 hsi: 9:00-12:00 / 13:00 - 16:00
2kospi: 9:00-12:00 /12:00 -15:00
3 mdax: 9:00-14:00/14:00-20:00
4 n100: 9:00-13:00/13:00-17:30
5aex: 9:00-13:00/13:00-17.30
6cac: 8:00-14:00/14:00-20:00  9:00 17:39
7dax: 9:00-14:00/14:00-20:00  3:00-11:30
8estx: 9:00-14:00/14:00-20:00 9:00-17:00
9xjo: 10:00-13:00/13:00~16:30
10nk200请从2010年开始，他的交易时间变了
nk200 9:00-11:00/12：30-15：00

添加了filling步骤

%}

clear

write_sel = false;
key_str = '主要指数测试';

dN = 'S40_america';
symbols = {'hsi','kospi','mdax','n100','aex','cac40','dax','estx50','xjo','nk200'};
t_cut = {[9,30,12,0;13,0,16,0],[9,30,12,0;12,0,15,0],[9,30,12,0;12,0,15,0],...
    [9,30,13,0;13,0,17,30],[9,30,13,0;13,0,17,30],[9,30,14,0;14,0,20,0],...
    [3,30,7,0;7,0,11,0],[9,30,12,0;12,0,17,0],[10,30,13,0;13,0,16,00],...
    [9,30,11,0;12,30,15,0]};
T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);

if write_sel
    obj_wd = wordcom(fullfile(pwd,sprintf('%s.doc',key_str)));
end

tradeDay_av = zeros(T_symbols,1);
re = cell(T_symbols,1);
for i_sym =  3%1:T_symbols
    sub_t_cut = t_cut{i_sym};
    sub_t1 = sub_t_cut(1,1:2);
    sub_t2 = sub_t_cut(1,3:4);
    sub_t3 = sub_t_cut(2,1:2);
    sub_t4 = sub_t_cut(2,3:4);
    
    sub_t0 = datenum(0,0,0,sub_t1(1),sub_t1(2),0);
    sub_t_dur = datenum(0,0,0,0,1,0);
    sub_t_t = datenum(0,0,0,sub_t2(1),sub_t2(2),0);
    sub_ind1 = sub_t0+sub_t_dur:sub_t_dur:sub_t_t;
    sub_ind1 = datevec(sub_ind1);
    sub_ind1 = sub_ind1(:,4)*100+sub_ind1(:,5);
    
    sub_t0 = datenum(0,0,0,sub_t3(1),sub_t3(2),0);
    sub_t_dur = datenum(0,0,0,0,1,0);
    sub_t_t = datenum(0,0,0,sub_t4(1),sub_t4(2),0);
    sub_ind2 = sub_t0+sub_t_dur:sub_t_dur:sub_t_t;
    sub_ind2 = datevec(sub_ind2);
    sub_ind2 = sub_ind2(:,4)*100+sub_ind2(:,5);
    
    sub_t_all = [sub_ind1;sub_ind2];
    T_all_t = length(sub_t_all);
    P = [];
    P.feeOpen=5/100000;
    P.feeClose=5/100000;
    P.matchRecord=1;%匹配数据源：沪深300
    P.tradeRecord=1;%交易数据源：股指期货主力合约
    P.tradeMin=length(sub_ind1);%使用早盘135分钟K线数据进行分形匹配
    P.dayMin=length(sub_t_all);%每个交易日共225根1分钟K线
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
    id_sel = (temp_t>sub_t1(1)*100+sub_t1(2)&temp_t<=sub_t2(1)*100+sub_t2(2)) | ...
        (temp_t>sub_t3(1)*100+sub_t3(2)&temp_t<=sub_t4(1)*100+sub_t4(2));
    x = x(id_sel,:);
    temp_t = temp_t(id_sel);
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
        if sum(sub_ind)<length(sub_t_all)-10 || ~eq(sub_t_all(1),sub_tt(1))
            ind_miss(sub_ind) = true;
            ind_miss_u(i) = true;
        else
           %filling 使用前面的值
            [~,ia,ib] = intersect(sub_t_all,sub_tt);
            sub_sub_x1 = nan(T_all_t,size(sub_x,2));
            sub_sub_x1(ia,:) = sub_x(ib,:);
            sub_sub_x1 = fillmissing(sub_sub_x1,'previous');
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