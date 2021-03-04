%{
1���۲ⴰ��������һ�������
2������β����close�۸�����
3��ͳ�Ƶڶ����������ȷ������
4���Եڶ���close�۸�����
5����������stop loss����
6�������ѻ���˫��ǧ��֮1.5

��Ҫ����޸Ŀ��
�ɶ�ɿ�
%}

clear

key_str = 'A�ɼ���-�Ľ�2';
symbols = fetchmysql('show tables from ycz_min_series',2);
T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
sta_re2 = sta_re;
Y1 = sta_re;
Y2 = sta_re;
error_ind = zeros(T_symbols,1);
%sql_str_r = 'select tradeDate,chgPct from yuqerdata.yq_dayprice where symbol = ''%s'' order by tradeDate';
sql_str_r = 'select tradeDate,closeprice/openprice-1 from yuqerdata.yq_dayprice where symbol = ''%s'' order by tradeDate';

parfor i_sym =1:T_symbols
    %��������
%     try 
        P = [];
        P.feeOpen=1.5/1000/2;
        P.feeClose=1.5/1000/2;
        P.matchRecord=1;%ƥ������Դ������300
        P.tradeRecord=1;%��������Դ����ָ�ڻ�������Լ
        P.tradeMin=240;%ʹ������120����K�����ݽ��з���ƥ��
        P.dayMin=240;%ÿ�������չ�240��1����K��
        P.M=20;%��M����Ϊ���ƵĽ�����
        P.muchPara=0.5;%�������ǻ��µ�����
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
        ind = SMTTradingModelTool_update1(closeprice,openprice,t,closeprice,openprice,t,P);
        r = fetchmysql(sprintf(sql_str_r,title_str(3:end)),2);                
        t1= cellstr(datestr(ind(:,1),'yyyy-mm-dd'));
        [sub_tref,ia,ib] = intersect(r(:,1),t1(:,1));
        x=[ind(ib,:),cell2mat(r(ia,end))];
        
        ind1 = x(:,2);
        ind1(2:end) = ind1(1:end-1);
        %ind1 = -ind1;
        ind1(ind1<0) = 0;
        r = x(:,end);
        r2 = r;
        fee_ind = find(~eq(diff(ind1),0))+1;
        r=r.*ind1;
        r2 = r2.*ind1;
        r2(fee_ind) = r2(fee_ind) -P.feeOpen-P.feeClose;
        y_c2 = cumprod(1+r);
        y_c = cumprod(1+r2);
        %ͳ�Ʋ���
        [v0,v_str0] = curve_static(y_c,[],false);
        [v,v_str] = ad_trans_sta_info(v0,v_str0); 
        result2 = [v_str;v]';
        result = [{'',title_str};result2];
        sta_re{i_sym} = result;
        sprintf('%s %d-%d',key_str,i_sym,T_symbols)
        Y1{i_sym} = [sub_tref,num2cell([y_c,y_c2])];
        %������
        ind1 = x(:,2);
        ind1(2:end) = ind1(1:end-1);
        ind1 = -ind1;
        ind1(ind1<0) = 0;
        r = x(:,end);
        r2 = r;
        r = r.*ind1;
        r2 = r2.*ind1;
        fee_ind = find(~eq(diff(ind1),0))+1;
        r2(fee_ind) = r2(fee_ind) -P.feeOpen-P.feeClose;
        y_c2 = cumprod(1+r);
        y_c = cumprod(1+r2);
        %ͳ�Ʋ���
        [v0,v_str0] = curve_static(y_c,[],false);
        [v,v_str] = ad_trans_sta_info(v0,v_str0); 
        result2 = [v_str;v]';
        result = [{'',title_str};result2];
        sta_re2{i_sym} = result;
        sprintf('%s %d-%d',key_str,i_sym,T_symbols)
        Y2{i_sym} = [sub_tref,num2cell([y_c,y_c2])];
        
        %     catch
%         error_ind(i_sym) = 1;
%         sprintf('Error %s %d-%d',key_str,i_sym,T_symbols)
%     end
%     
end
y = [sta_re{:}];
y = y(:,[1,2:2:end]);
y = y';

y2 = [sta_re2{:}];
y2 = y2(:,[1,2:2:end]);
y2 = y2';


