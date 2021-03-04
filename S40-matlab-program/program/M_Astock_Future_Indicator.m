clear

key_str = 'S40SMT���Թ�ָ�ڻ�����';
symbols = {'ZJIH','ZJIC','ZJIF'};
symbols_tdx = {'ccfx_ih','ccfx_ic','ccfx_if'};
dN_tdx = 'Future_min_data';
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
for i_sym = 1:T_symbols
    %��������
    P = [];
    P.feeOpen=5/100000;
    P.feeClose=5/100000;
    P.matchRecord=1;%ƥ������Դ������300
    P.tradeRecord=1;%��������Դ����ָ�ڻ�������Լ
    P.tradeMin=120;%ʹ������120����K�����ݽ��з���ƥ��
    P.dayMin=240;%ÿ�������չ�240��1����K��
    P.M=20;%��M����Ϊ���ƵĽ�����
    P.muchPara=0.5;%�������ǻ��µ�����
    P.deanMethod=3;%1���ϵ��/2ŷʽ����/3���Ͼ���/4�����پ���
    P.stopMethod=3;%1���̼�ֹ��/2������ֹ�𣬲������ܷ���/3�������̼�ֹ�𣬷��򴥷���ֹ��
    P.testStart='2010-4-16';
    P.trade_mode = 2;%1ֻ��� 2 ��ֺͿղ�

    title_str = symbols{i_sym};
    sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
        'hour(tradingdate),minute(tradingdate),open,close from S40.%s'];
    sub_sql_str = sprintf(sql_str,title_str);
    x = fetchmysql(sub_sql_str);
    if isempty(x)
        sql_str_tdx = ['select t_year,t_month,t_day,',...
            't_hour,t_minute,open,close from %s.%s ',...
            ' order by tradingdate'];
        sub_sql_str = sprintf(sql_str_tdx,dN_tdx,symbols_tdx{i_sym});
        x = fetchmysql(sub_sql_str);
    else
        sub_t_max = x(end,1:6);
        sub_t_max(end) = 0;
        sub_t_max = datestr(datenum(sub_t_max),'yyyy-mm-dd HH:MM:SS');
        sql_str_tdx = ['select t_year,t_month,t_day,',...
            't_hour,t_minute,open,close from %s.%s ',...
            'where tradingdate>''%s'' order by tradingdate'];
        sub_sql_str = sprintf(sql_str_tdx,dN_tdx,symbols_tdx{i_sym},sub_t_max);
        x1 = fetchmysql(sub_sql_str);
        x = cat(1,x,x1);
    end
    if isempty(x)
        sprintf('%s %s ʱ�䲻��6�꣬����',key_str,title_str)
        continue
    end
    temp = x(:,4)*100+x(:,5);
    id_sel = temp>930&temp<=1500;
    x = x(id_sel,:);

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
    min_day_num = 210*4;
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
    ah=gca;
    title(ah,title_str)
    if write_sel
        obj_wd.pasteFigure(h,title_str);
    end
    y_c = cumprod(tradeYield(:,2)+1);
    %ͳ�Ʋ���
    [v0,v_str0] = curve_static(y_c,[],false);
    [v,v_str] = ad_trans_sta_info(v0,v_str0); 
    result2 = [v_str;v]';
    result = [{'',title_str};[result1;result2]];
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