clear

key_str = '指数测试';

dns = {'IF','IH','IC';'300','50','500'};

T_symbols = size(dns,2);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);
for i_sym = 1%1:T_symbols
    %参数设置
    %try 
        P = [];
        P.feeOpen=5/100000;
        P.feeClose=5/100000;
        P.matchRecord=1;%匹配数据源：沪深300
        P.tradeRecord=2;%交易数据源：股指期货主力合约
        P.tradeMin=120;%使用早盘120分钟K线数据进行分形匹配
        P.dayMin=240;%每个交易日共240根1分钟K线
        P.M=20;%找M个最为相似的交易日
        P.muchPara=0.5;%多数上涨或下跌比例
        P.deanMethod=3;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离
        P.stopMethod=3;%1收盘价止损/2触发价止损，不考虑能否交易/3跳开则开盘价止损，否则触发价止损
        P.testStart='2010-4-16';

        title_str = dns{1,i_sym};
        sub_tn1 = dns{1,i_sym};
        sub_tn2 = dns{2,i_sym};
        
        %期货数据
        sql_str1 = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
            'hour(tradingdate),minute(tradingdate),openprice,closeprice from  S28.wind_%s order by tradingdate'];
        sub_x1 = fetchmysql(sprintf(sql_str1,sub_tn1));        
        %指数数据
        sql_str2 = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
            'hour(tradingdate),minute(tradingdate),openprice,closeprice from  S28.wind_%s order by tradingdate'];
        sub_x2 = fetchmysql(sprintf(sql_str2,sub_tn2));     
        
        [t1,testStart1,openprice1,closeprice1] = S40_preprocessingdata(sub_x1);
        [t2,testStart2,openprice2,closeprice2] = S40_preprocessingdata(sub_x2);
        
        P.testStart = testStart1;
        
        %SMTTradingModelTool(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate)
        [tradeYield,result1,tradeDetail,yearDetail,h] = SMTTradingModelTool(closeprice2,openprice2,t2,closeprice1,openprice1,t1,P);
        y_c = cumprod(tradeYield(:,2)+1);
        %统计参数
        [v0,v_str0] = curve_static(y_c,[],false);
        [v,v_str] = ad_trans_sta_info(v0,v_str0); 
        result2 = [v_str;v]';
        result = [{'',title_str};[result1;result2]];
        sta_re{i_sym} = result;
%     catch
%         error_ind(i_sym) = 1;
%     end
    sprintf('%s %d-%d',key_str,i_sym,T_symbols)
end
y = [sta_re{:}];
y = y(:,[1,2:2:end]);
y = y';
