%{
��Ʒ�ڻ�����
�ѵ�
��ҹ�̽��ף���δ���
��������ͬ����δ��������ж��Ƿ���ֹ�𿿵��ǵ�����
��ʱ������ ҹ�̲��迼��
�ж��Ƿ�ֹ���������

����������
1������Դ
2��ͨ����

���ڻ�-������Ϊʮ���֮�壬���Ϊʮ���֮1
%}

clear

info = ['select secFullName,contractObject,exchangeCD  ',...
    'from yuqerdata.yq_FutuGet group by contractObject order by contractObject'];
info = fetchmysql(info,2);

key_str = 'S40SMT���Թ�����Ʒ�ڻز���';

dN = 'Future_min_Jinshuyuan';
dN_tdx = 'future_min_data';
symbols = fetchmysql(sprintf('show tables from %s',dN),2);
T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
error_ind = zeros(T_symbols,1);

write_sel = true;
if write_sel
    pn_write = fullfile(pwd,'������');
    if ~exist(pn_write,'dir')
        mkdir(pn_write)
    end
    obj_wd = wordcom(fullfile(pn_write,sprintf('%s.doc',key_str)));
    xls_fn = fullfile(pn_write,sprintf('%s.xlsx',key_str));
end

tradeDay_av = zeros(T_symbols,1);
for i_sym = 1:T_symbols

    P = [];
    P.feeOpen=5/100000;
    P.feeClose=5/100000;
    P.matchRecord=1;%ƥ������Դ������300
    P.tradeRecord=1;%��������Դ����ָ�ڻ�������Լ
    P.tradeMin=135;%ʹ������135����K�����ݽ��з���ƥ��
    P.dayMin=225;%ÿ�������չ�225��1����K��
    P.M=20;%��M����Ϊ���ƵĽ�����
    P.muchPara=0.5;%�������ǻ��µ�����
    P.deanMethod=3;%1���ϵ��/2ŷʽ����/3���Ͼ���/4�����پ���
    P.stopMethod=3;%1���̼�ֹ��/2������ֹ�𣬲������ܷ���/3�������̼�ֹ�𣬷��򴥷���ֹ��
    P.testStart='2010-4-16';
    P.trade_mode = 2;%1ֻ��� 2 ��ֺͿղ�

    title_str = symbols{i_sym};
    temp_id = strcmpi(info(:,2),title_str);
    if any(temp_id)
        title_str2 = info{temp_id,1};
        tn_tdx = sprintf('%s_%s',info{temp_id,3},info{temp_id,2});
    else
        title_str2 = title_str;
        tn_tdx = [];
    end
    
    sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
        'hour(tradingdate),minute(tradingdate),openprice,closeprice from %s.%s ',...
        'order by tradingdate'];
    sub_sql_str = sprintf(sql_str,dN,title_str);
    x = fetchmysql(sub_sql_str);
    if ~isempty(tn_tdx)
        sub_t_max = x(end,1:6);
        sub_t_max(end) = 0;
        sub_t_max = datestr(datenum(sub_t_max),'yyyy-mm-dd HH:MM:SS');
        sql_str_tdx = ['select t_year,t_month,t_day,',...
            't_hour,t_minute,open,close from %s.%s ',...
            'where tradingdate>''%s'' order by tradingdate'];
        sub_sql_str = sprintf(sql_str_tdx,dN_tdx,tn_tdx,sub_t_max);
        x1 = fetchmysql(sub_sql_str);
        x = cat(1,x,x1);
    else
        sql_str_tdx = ['select t_year,t_month,t_day,',...
            't_hour,t_minute,open,close from %s.%s ',...
            ' order by tradingdate'];
        sub_sql_str = sprintf(sql_str_tdx,dN_tdx,tn_tdx);
        x = fetchmysql(sub_sql_str);
    end
    if isempty(x)
        sprintf('%s %s ʱ�䲻��6�꣬����',key_str,title_str)
        continue
    end
    temp = x(:,4)*100+x(:,5);
    id_sel = temp>900&temp<=1500;
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
    [tradeYield,result1,tradeDetail,yearDetail,h] = SMTTradingModelTool2(closeprice,openprice,t,closeprice,openprice,t,P);
    ah = gca;
    title(ah,title_str2);
    if write_sel
        obj_wd.pasteFigure(h,title_str2);
    end
    y_c = cumprod(tradeYield(:,2)+1);
    %ͳ�Ʋ���
    [v0,v_str0] = curve_static(y_c,[],false);
    [v,v_str] = ad_trans_sta_info(v0,v_str0); 
    result2 = [v_str;v]';
    result = [{'',title_str2};[result1;result2]];
    sta_re{i_sym} = result;
    sprintf('%s %d-%d',key_str,i_sym,T_symbols)
    
end
y = [sta_re{:}];
y = y(:,[1,2:2:end]);
y = y';
if write_sel
    obj_wd.CloseWord()
    xlswrite(xls_fn,y);
end
