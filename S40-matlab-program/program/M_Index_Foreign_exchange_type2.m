%{
����
���޼����������ʣ���
audcad,audchf,audjpy,audnzd,audusd,nzdusd,nzdjpy,euraud,eurnzd,usdjpy
8:00am to 14:00 �Ǽ������ʱ���
14.00am to 20:00 �ǽ���ʱ���
%}

clear

write_sel = false;
key_str = '�������';

dN = 'foreign_index_min_V2';
symbols = fetchmysql(sprintf('show tables from %s',dN),2);
symbols1 = strsplit('audcad,audchf,audjpy,audnzd,audusd,nzdusd,nzdjpy,euraud,eurnzd,usdjpy',',');
symbols = setdiff(symbols,symbols1);
T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);

if write_sel
    obj_wd = wordcom(fullfile(pwd,sprintf('%s.doc',key_str)));
end

tradeDay_av = zeros(T_symbols,1);
re = cell(T_symbols,1);

parfor i_sym =  1:T_symbols

    P = [];
    P.feeOpen=1/100000;
    P.feeClose=1/100000;
    P.matchRecord=1;%ƥ������Դ������300
    P.tradeRecord=1;%��������Դ����ָ�ڻ�������Լ
    P.tradeMin=360;%ʹ������135����K�����ݽ��з���ƥ��
    P.dayMin=360;%ÿ�������չ�225��1����K��
    P.M=20;%��M����Ϊ���ƵĽ�����
    P.muchPara=0.5;%�������ǻ��µ�����
    P.deanMethod=3;%1���ϵ��/2ŷʽ����/3���Ͼ���/4�����پ���
    P.stopMethod=3;%1���̼�ֹ��/2������ֹ�𣬲������ܷ���/3�������̼�ֹ�𣬷��򴥷���ֹ��
    P.testStart='2010-4-16';
    P.trade_mode = 2;%1ֻ��� 2 ��ֺͿղ�

    title_str = symbols{i_sym};  
    sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
        'hour(tradingdate),minute(tradingdate),open,close from %s.%s ',...
        ' order by tradingdate'];
    sub_sql_str = sprintf(sql_str,dN,title_str);
    x = fetchmysql(sub_sql_str);
    if isempty(x)
        sprintf('%s %s ʱ�䲻��4�꣬����',key_str,title_str)
        continue
    end
    temp = x(:,4)*100+x(:,5);
    id_sel = temp>800&temp<=2000 ;
    %id_sel = (temp>=800&temp<=1000) |(temp>=1100&temp<=1430);
    x = x(id_sel,:);
    if isempty(x)
        continue
    end
    t = datenum([x(:,1:5),zeros(size(x(:,1)))]);
    %ͳ���м���ͣ�Ƶ���������޳�
    day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
    day_tick_u = unique(day_tick);
    ind_miss = false(size(day_tick));
    ind_miss_u = false(size(day_tick_u));
    T = length(day_tick_u);
    temp_tradedays = zeros(T,1);
    for i = 1:T
        sub_ind = eq(day_tick,day_tick_u(i));
        if ~eq(sum(sub_ind),720)
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

    %��ʼʱ���趨
    min_day_num = 210*4;
    if length(day_tick_u)<min_day_num
        sprintf('%s %s ʱ�䲻��4�꣬����',key_str,title_str)
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
%     if write_sel
%         obj_wd.pasteFigure(h,title_str);
%     end
    y_c = cumprod(tradeYield(:,2)+1);
    %ͳ�Ʋ���
    [v0,v_str0] = curve_static(y_c,[],false);
    [v,v_str] = ad_trans_sta_info(v0,v_str0); 
    result2 = [v_str;v]';
    result = [{'',title_str};[result1;result2]];
    sta_re{i_sym} = result;
    sprintf('%s %d-%d',key_str,i_sym,T_symbols)
    re{i_sym} = tradeYield;
end
if write_sel
    obj_wd.CloseWord()
end
y = [sta_re{:}];
y = y(:,[1,2:2:end]);
y = y';
xlswrite(sprintf('%s.xlsx',key_str),y);

save('foreign_exchange1_u2.mat','re','symbols')