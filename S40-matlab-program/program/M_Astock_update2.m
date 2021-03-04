%{
��������Ľ�һ��
��Ʊ-ֻ����
1������������������ݼ������ƶȵ�ʱ�򣬵õ�ǰn�������Ƶ����ӣ������������棬ֻ�е���ƽ�������棨������С���棩����һ����ֵ������>ǧ��֮1.5����ʱ����ܴ����źš�
2�����ھ�������¶��壺ԭ���߶�����3�־��롣���Զ���DTW��dynamic time wrapper)Ϊ�����־��롣
���������Ϥ����������һ�¡������ʵ�ֿ�������github��������һ��matlab�����кܶ�python�ġ�
https://github.com/xdjcl/DTW

%}
close all
clear

key_str = 'A�ɼ���-�Ľ�2';
symbols = fetchmysql('show tables from ycz_min_series',2);
T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);
for i_sym =1% 1:T_symbols
    %��������
    %try 
        P = [];
        P.feeOpen=1.5/1000/2;
        P.feeClose=1.5/1000/2;
        P.matchRecord=1;%ƥ������Դ������300
        P.tradeRecord=1;%��������Դ����ָ�ڻ�������Լ
        P.tradeMin=120;%ʹ������120����K�����ݽ��з���ƥ��
        P.dayMin=240;%ÿ�������չ�240��1����K��
        P.M=20;%��M����Ϊ���ƵĽ�����
        P.muchPara=0.7;%�������ǻ��µ�����
        P.deanMethod=3;%1���ϵ��/2ŷʽ����/3���Ͼ���/4�����پ���/5 dynamic time warping
        P.stopMethod=3;%1���̼�ֹ��/2������ֹ�𣬲������ܷ���/3�������̼�ֹ�𣬷��򴥷���ֹ��
        P.testStart='2010-4-16';
        P.trade_mode = 1;%1ֻ��� 2 ��ֺͿղ�
        P.cut_return = -inf;
        
        title_str = symbols{i_sym};
        sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
            'hour(tradingdate),minute(tradingdate),open,close from ycz_min_series.%s'];
        sub_sql_str = sprintf(sql_str,title_str);
        x = fetchmysql(sub_sql_str);
        if isempty(x)
            sprintf('%s %s ʱ�䲻��6�꣬����',key_str,title_str)
            continue
        end
        t = datenum([x(:,1:5),zeros(size(x(:,1)))]);
        %ͳ���м���ͣ�Ƶ���������޳�
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

        %��ʼʱ���趨
        min_day_num = 210*6;
        if length(day_tick_u)<min_day_num
            sprintf('%s %s ʱ�䲻��6�꣬����',key_str,title_str)
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
        %ͳ�Ʋ���
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
