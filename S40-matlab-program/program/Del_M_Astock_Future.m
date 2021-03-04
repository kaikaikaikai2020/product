clear

key_str = 'ָ������';

dns = {'IF','IH','IC';'300','50','500'};

T_symbols = size(dns,2);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);
for i_sym = 1%1:T_symbols
    %��������
    %try 
        P = [];
        P.feeOpen=5/100000;
        P.feeClose=5/100000;
        P.matchRecord=1;%ƥ������Դ������300
        P.tradeRecord=2;%��������Դ����ָ�ڻ�������Լ
        P.tradeMin=120;%ʹ������120����K�����ݽ��з���ƥ��
        P.dayMin=240;%ÿ�������չ�240��1����K��
        P.M=20;%��M����Ϊ���ƵĽ�����
        P.muchPara=0.5;%�������ǻ��µ�����
        P.deanMethod=3;%1���ϵ��/2ŷʽ����/3���Ͼ���/4�����پ���
        P.stopMethod=3;%1���̼�ֹ��/2������ֹ�𣬲������ܷ���/3�������̼�ֹ�𣬷��򴥷���ֹ��
        P.testStart='2010-4-16';

        title_str = dns{1,i_sym};
        sub_tn1 = dns{1,i_sym};
        sub_tn2 = dns{2,i_sym};
        
        %�ڻ�����
        sql_str1 = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
            'hour(tradingdate),minute(tradingdate),openprice,closeprice from  S28.wind_%s order by tradingdate'];
        sub_x1 = fetchmysql(sprintf(sql_str1,sub_tn1));        
        %ָ������
        sql_str2 = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
            'hour(tradingdate),minute(tradingdate),openprice,closeprice from  S28.wind_%s order by tradingdate'];
        sub_x2 = fetchmysql(sprintf(sql_str2,sub_tn2));     
        
        [t1,testStart1,openprice1,closeprice1] = S40_preprocessingdata(sub_x1);
        [t2,testStart2,openprice2,closeprice2] = S40_preprocessingdata(sub_x2);
        
        P.testStart = testStart1;
        
        %SMTTradingModelTool(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate)
        [tradeYield,result1,tradeDetail,yearDetail,h] = SMTTradingModelTool(closeprice2,openprice2,t2,closeprice1,openprice1,t1,P);
        y_c = cumprod(tradeYield(:,2)+1);
        %ͳ�Ʋ���
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
