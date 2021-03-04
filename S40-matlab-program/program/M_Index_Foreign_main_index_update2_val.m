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
升级
a）找出所有交易日的开始时间的中位值（绝大多数都是这个时间点开始有数据）A
b）找出所有交易日收盘时间的中位值（绝大多数都是这个时间点没有数据）B
c）过滤时间点不在A-B的日期
d）A-B中间如果有少量数据缺失（10分钟以内），使用前一时间filling；
e）以A-B中间时间为分割点，前半部分计算相关系数，后半部分执行；

画出指数，验证

%}

clear

write_sel = false;
key_str = '国外主要指数-验证';

dN = 'S40_america';
symbols = {'hsi','kospi','mdax','n100','aex','cac40','dax','estx50','xjo','nk200'};
t_cut = {[9,30,12,0;13,0,16,0],[9,30,12,0;12,0,15,0],[9,30,12,0;12,0,15,0],...
    [9,30,13,0;13,0,17,30],[9,30,13,0;13,0,17,30],[9,30,14,0;14,0,20,0],...
    [3,30,7,0;7,0,11,0],[9,30,12,0;12,0,17,0],[10,30,13,0;13,0,16,00],...
    [9,30,11,0;12,30,15,0]};
T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);


sub_t0 = datenum(0,0,0,0,0,0);
sub_t_dur = datenum(0,0,0,0,1,0);
sub_t_t = datenum(0,0,0,23,59,0);
sub_ind1 = sub_t0+sub_t_dur:sub_t_dur:sub_t_t;
sub_ind1 = datevec(sub_ind1);
day_minute = sub_ind1(:,4)*100+sub_ind1(:,5);

if write_sel
    obj_wd = wordcom(fullfile(pwd,sprintf('%s.doc',key_str)));
end

tradeDay_av = zeros(T_symbols,1);
re = cell(T_symbols,1);
re2 = re;
for i_sym =  1:T_symbols
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
    %这里增加判断时间
    day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
    day_tick_u = unique(day_tick);
    t_min = zeros(size(day_tick_u));
    t_max = t_min;
    for i = 1:length(t_min)
        sub_ind = eq(day_tick,day_tick_u(i));
        sub_t = temp_t(sub_ind);
        t_min(i) = min(sub_t);
        t_max(i) = max(sub_t);        
    end
    t_min1 = median(t_min);
    t_max1 = median(t_max);
    
    sub_t_all = day_minute(day_minute>=t_min1 & day_minute<=t_max1);  
    T_all_t = length(sub_t_all);
    id_sel = temp_t>=t_min1&temp_t<=t_max1;
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
    re2{i_sym} = [length(temp_tradedays),t_min1,t_max1]';
    figure
    plot(x(:,end));
    t_str = cellstr(datestr(t,'yyyymmdd'));
    T = length(t_str);
    set(gca,'xlim',[0,T]);
    set(gca,'XTick',floor(linspace(1,T,10)));
    set(gca,'XTickLabel',t_str(floor(linspace(1,T,10))));
    set(gca,'XTickLabelRotation',90)
    title(title_str);
%     openprice = x(:,end-1);
%     closeprice = x(:,end);
%     %SMTTradingModelTool(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate)
%     P.tradeMin=floor(length(sub_t_all)/2);%使用早盘135分钟K线数据进行分形匹配
%     P.dayMin=length(sub_t_all);%每个交易日共225根1分钟K线
%     [tradeYield,result1,tradeDetail,yearDetail,h,assure_ratio] = SMTTradingModelTool(closeprice,openprice,t,closeprice,openprice,t,P);
%     re{i_sym} = [tradeYield,assure_ratio];
%     ah = gca;
%     title(ah,title_str);
%     if write_sel
%         obj_wd.pasteFigure(h,title_str);
%     end
%     y_c = cumprod(tradeYield(:,2)+1);
%     %统计参数
%     [v0,v_str0] = curve_static(y_c,[],false);
%     [v,v_str] = ad_trans_sta_info(v0,v_str0); 
%     result2 = [v_str;v]';
%     result = [{'',title_str};[result1;result2]];
%     sta_re{i_sym} = result;
%     sprintf('%s %d-%d',key_str,i_sym,T_symbols)
    
end
if write_sel
    obj_wd.CloseWord()
end
y = [sta_re{:}];
y = y(:,[1,2:2:end]);
y = y';
xlswrite(sprintf('%s.xlsx',key_str),y);

save('foreign_main_index.mat','re','symbols')