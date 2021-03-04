%{
两个方面改进一下
股票-只做多
1）当我们用上午的数据计算相似度的时候，得到前n个最相似的日子，计算下午收益，只有当其平均的收益（或者最小收益）超过一定阀值（至少>千分之1.5）的时候才能触发信号。
2）对于距离的重新定义：原作者定义了3种距离。可以定义DTW（dynamic time wrapper)为第四种距离。
如果您不熟悉，可以搜索一下。具体的实现可以搜索github，下面有一个matlab。还有很多python的。
https://github.com/xdjcl/DTW

%}
close all
clear

key_str = 'A股计算-改进2';
symbols = fetchmysql('show tables from ycz_min_series',2);
T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);
for i_sym =1% 1:T_symbols
    %参数设置
    %try 
        P = [];
        P.feeOpen=1.5/1000/2;
        P.feeClose=1.5/1000/2;
        P.matchRecord=1;%匹配数据源：沪深300
        P.tradeRecord=1;%交易数据源：股指期货主力合约
        P.tradeMin=120;%使用早盘120分钟K线数据进行分形匹配
        P.dayMin=240;%每个交易日共240根1分钟K线
        P.M=20;%找M个最为相似的交易日
        P.muchPara=0.7;%多数上涨或下跌比例
        P.deanMethod=3;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离/5 dynamic time warping
        P.stopMethod=3;%1收盘价止损/2触发价止损，不考虑能否交易/3跳开则开盘价止损，否则触发价止损
        P.testStart='2010-4-16';
        P.trade_mode = 1;%1只多仓 2 多仓和空仓
        P.cut_return = -inf;
        
        title_str = symbols{i_sym};
        sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
            'hour(tradingdate),minute(tradingdate),open,close from ycz_min_series.%s'];
        sub_sql_str = sprintf(sql_str,title_str);
        x = fetchmysql(sub_sql_str);
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
        [tradeYield,result1,tradeDetail,yearDetail,h] = SMTTradingModelTool_update2(closeprice,openprice,t,closeprice,openprice,t,P);
        y_c = cumprod(tradeYield(:,2)+1);
        %统计参数
        [v0,v_str0] = curve_static(y_c,[],false);
        [v,v_str] = ad_trans_sta_info(v0,v_str0); 
        result2 = [v_str;v]';
        result = [{'',title_str};[result1;result2]];
        sta_re{i_sym} = result;
        sprintf('%s %d-%d',key_str,i_sym,T_symbols)
%     catch
%         error_ind(i_sym) = 1;
%         sprintf('Error %s %d-%d',key_str,i_sym,T_symbols)
%     end
    
end
% y = [sta_re{:}];
% y = y(:,[1,2:2:end]);
% y = y';
