clear
%{
每次算出来的概率作为权重，比如计算出来0.6的可能性涨，我们就以0.6作为权重。。。
如果下一个股票0.55的概率涨，权重就为0.55. 升级版本
%}
key_str = '指数测试';
symbols = {'sh000016','sh000905','sh000300'};
f_symbols = {'ZJIH','ZJIC','ZJIF'};
T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);
for i_sym = 1:T_symbols
    %参数设置
    try 
        P = [];
        P.feeOpen=5/100000;
        P.feeClose=5/100000;
        P.matchRecord=1;%匹配数据源：沪深300
        P.tradeRecord=1;%交易数据源：股指期货主力合约
        P.tradeMin=120;%使用早盘120分钟K线数据进行分形匹配
        P.dayMin=240;%每个交易日共240根1分钟K线
        P.M=20;%找M个最为相似的交易日
        P.muchPara=0.5;%多数上涨或下跌比例
        P.deanMethod=3;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离
        P.stopMethod=3;%1收盘价止损/2触发价止损，不考虑能否交易/3跳开则开盘价止损，否则触发价止损
        P.testStart='2010-4-16';
        P.trade_mode = 2;%1只多仓 2 多仓和空仓

        title_str = symbols{i_sym};
        sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
            'hour(tradingdate),minute(tradingdate),open,close from S40.%s'];
        sub_sql_str = sprintf(sql_str,title_str);
        x = fetchmysql(sub_sql_str);
        temp = x(:,4)*100+x(:,5);
        id_sel = temp>930&temp<=1500;
        x = x(id_sel,:);
        if isempty(x)
            sprintf('%s %s 时间不足6年，跳出',key_str,title_str)
            continue
        end
        t = datenum([x(:,1:5),zeros(size(x(:,1)))]);
        %统计中间有停牌的情况，并剔除
        day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
        day_tick_u = unique(day_tick);
        ind_miss = false(size(day_tick));
        ind_miss_u = false(size(day_tick_u));
        T = length(day_tick_u);
        for i = 1:T
            sub_ind = eq(day_tick,day_tick_u(i));
            if sum(sub_ind)<240
                ind_miss(sub_ind) = true;
                ind_miss_u(i) = true;
            end
        end

        day_tick_u(ind_miss_u) = [];
        x(ind_miss,:) = [];
        t(ind_miss,:) = [];

        %初始时间设定
        min_day_num = 210*6;
        if length(day_tick_u)<min_day_num
            sprintf('%s %s 时间不足6年，跳出',key_str,title_str)
            continue
        else
            temp = num2str(day_tick_u(min_day_num/2+1));
            P.testStart=sprintf('%s-%s-%s',temp(1:4),temp(5:6),temp(7:8));
        end

        openprice = x(:,end-1);
        closeprice = x(:,end);
        %SMTTradingModelTool(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate)
        [tradeYield,result1,tradeDetail,yearDetail,h] = SMTTradingModelTool_update_weight(closeprice,openprice,t,closeprice,openprice,t,P);
        ah=gca;
        title(ah,title_str)
        y_c = cumprod(tradeYield(:,2)+1);
        %统计参数
        [v0,v_str0] = curve_static(y_c,[],false);
        [v,v_str] = ad_trans_sta_info(v0,v_str0); 
        result2 = [v_str;v]';
        result = [{'',title_str};[result1;result2]];
        sta_re{i_sym} = result;
        sprintf('%s %d-%d',key_str,i_sym,T_symbols)
    catch
        error_ind(i_sym) = 1;
        sprintf('Error %s %d-%d',key_str,i_sym,T_symbols)
    end
    
end
y = [sta_re{:}];
y = y(:,[1,2:2:end]);
y = y';
