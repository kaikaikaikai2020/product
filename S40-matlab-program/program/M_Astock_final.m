%{
1���۲ⴰ��������һ�������
2������β����close�۸�����
3��ͳ�Ƶڶ����������ȷ������
4���Եڶ���close�۸�����
5����������stop loss����
6�������ѻ���˫��ǧ��֮1.5

�Խ�д�����ݿ�
%}

clear

tn = 'S40.A_stock_signal';
var_info={'symbol', 'tradenum', 'f_val'};

key_str = 'S40SMT����A�ɼ���';
symbols = fetchmysql('show tables from ycz_min_series',2);
T_symbols = length(symbols);
sta_re = cell(T_symbols,1);
sta_re2 = sta_re;
Y1 = sta_re;
Y2 = sta_re;
error_ind = zeros(T_symbols,1);
%sql_str_r = 'select tradeDate,chgPct from yuqerdata.yq_dayprice where symbol = ''%s'' order by tradeDate';
sql_str_r = 'select tradeDate,closeprice/openprice-1 from yuqerdata.yq_dayprice where symbol = ''%s'' order by tradeDate';
X = cell(T_symbols,1);
for i_sym =1:T_symbols
    %��������
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
        
        sub_sql_str =[ 'select tradenum,f_val from S40.A_stock_signal where ',...
            'symbol = ''%s'' order by tradenum'];
        ind_before = fetchmysql(sprintf(sub_sql_str,symbols{i_sym}));
        if isempty(ind_before)        
            P.testStart=sprintf('%s-%s-%s',temp(1:4),temp(5:6),temp(7:8));
        else
            P.testStart=datestr(ind_before(end,1)+1,'yyyy-mm-dd');
        end
    end

    openprice = x(:,end-1);
    closeprice = x(:,end);
    %SMTTradingModelTool(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate)
    if datenum(P.testStart)<=t(end)
        ind0 = SMTTradingModelTool_update1(closeprice,openprice,t,closeprice,openprice,t,P);
    else
        ind0 = [];
    end
    ind = [ind_before;ind0];
    %mysql_data
    if ~isempty(ind0)
        ind0 = num2cell(ind0(:,[1,1:end]));
        ind0(:,1) = symbols(i_sym);
        X{i_sym} = ind0';
    end
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
    r2(fee_ind) = r2(fee_ind) -P.feeOpen-P.feeClose;
    y_c2 = cumprod(1+r.*ind1);
    y_c = cumprod(1+r2.*ind1);
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
    fee_ind = find(~eq(diff(ind1),0))+1;
    r2(fee_ind) = r2(fee_ind) -P.feeOpen-P.feeClose;
    y_c2 = cumprod(1+r.*ind1);
    y_c = cumprod(1+r2.*ind1);
    %ͳ�Ʋ���
    [v0,v_str0] = curve_static(y_c,[],false);
    [v,v_str] = ad_trans_sta_info(v0,v_str0); 
    result2 = [v_str;v]';
    result = [{'',title_str};result2];
    sta_re2{i_sym} = result;
    sprintf('%s %d-%d',key_str,i_sym,T_symbols)
    Y2{i_sym} = [sub_tref,num2cell([y_c,y_c2])];          
end

X1 = [X{:}]';
datainsert_adair(tn,var_info,X1);

y = [sta_re{:}];
y = y(:,[1,2:2:end]);
y = y';

y2 = [sta_re2{:}];
y2 = y2(:,[1,2:2:end]);
y2 = y2';

write_sel = true;
if write_sel
    pn_write = fullfile(pwd,'������');
    if ~exist(pn_write,'dir')
        mkdir(pn_write)
    end
    obj_wd = wordcom(fullfile(pn_write,sprintf('%s.doc',key_str)));
    xls_fn = fullfile(pn_write,sprintf('%s.xlsx',key_str));
end

T = length(Y1);
tref = yq_methods.get_tradingdate();
Y = zeros(length(tref),T);
Y_0 = Y;
for i = 1:T
    sub_x = Y1{i};
    if ~isempty(sub_x)
        [~,ia,ib] = intersect(tref,sub_x(:,1));
        temp = cell2mat(sub_x(:,2:3));
        temp(2:end,:) = temp(2:end,:)./temp(1:end-1,:)-1;
        temp(1,:) = 0;
        Y(ia,i) = temp(ib,1);
        Y_0(ia,i) = temp(ib,2);
    end        
    sprintf('%d-%d',i,T)
end

sel_ind1 = sum(abs(Y),1)>1;
symbols2 = symbols(sel_ind1);
Y = Y(:,sel_ind1);
Y_0 = Y_0(:,sel_ind1);
sel_ind2 = sum(abs(Y),2)>1;
tref2 = tref(sel_ind2);
Y = Y(sel_ind2,:);
Y_0 = Y_0(sel_ind2,:);

symbols2 = cellfun(@(x) x(3:end),symbols2,'UniformOutput',false);

t2 = datenum(tref2);
index_pool = {'000300','000905','000016','000001'};
index_name = {'����300','��֤500','��֤50','��֤��ָ'};
sql_str1 = 'select tradeDate,closeIndex from   yuqerdata.yq_index where symbol = ''%s'' order by tradeDate';
sta_re = cell(size(index_name));
for index_sel = 1:length(index_pool)
    title_str = index_name{index_sel};
    sub_x = fetchmysql(sprintf(sql_str1,index_pool{index_sel}),2);
    sub_symbols = yq_methods.get_index_pool(index_pool{index_sel},sub_x{end,1});

    %���ƹ�Ʊ��
    [~,ia] = intersect(symbols2,sub_symbols);
    [sub_tref,ia1,ib1] = intersect(tref2,sub_x(:,1));

    sub_Y = Y(ia1,ia);
    sub_Y0 = Y_0(ia1,ia);
    sub_t_num = t2(ia1);
    sub_x = sub_x(ib1,:);
    sub_x_c = cell2mat(sub_x(:,2));
    sub_x_c = sub_x_c./sub_x_c(1);

    t_str = sub_tref;
    T = length(t_str);

    r  = mean(sub_Y,2);
    r0 = mean(sub_Y0,2);
    r(1) = 0;
    r0(1) = 0;
    r_c = cumprod(1+r);
    r_c0 = cumprod(1+r0);

    h1=figure;
    subplot(2,1,1)
    plot([r_c,r_c0,sub_x_c],'-','LineWidth',2);
    set(gca,'xlim',[0,T]);
    set(gca,'XTick',floor(linspace(1,T,15)));
    set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
    set(gca,'XTickLabelRotation',90)    
    setpixelposition(h1,[223,365,1345,600]);
    legend({'���������������','���������������',index_pool{index_sel}},'NumColumns',3,'Location','best')
    box off
    title(sprintf('%s-�������',title_str))
    %ÿ������
    t_y = year(sub_t_num);
    t_y_u = unique(t_y);
    r_year = zeros(size(t_y_u));
    for j = 1:length(t_y_u)
    sub_r = r(eq(t_y,t_y_u(j)),:);
    temp = cumprod(1+mean(sub_r,2));
    r_year(j) = temp(end)-1;
    end
    subplot(2,1,2)
    bar(t_y_u,r_year)
    box off
    title(sprintf('%s-ÿ������ͳ��',title_str))
    if write_sel
        obj_wd.pasteFigure(h1,title_str);
    end
    %ͳ�Ʋ���
    sub_re = cell(2,1);
    sub_y = [sub_x_c,r_c];
    sub_title_str = {title_str,sprintf('%s���',title_str)};
    for j = 1:2
        [v0,v_str0] = curve_static(sub_y(:,j),[],false);
        [v,v_str] = ad_trans_sta_info(v0,v_str0);
        if eq(j,1)
            sub_re{j} = [[{''},sub_title_str(j)];[v_str;v]'];
        else
            sub_re{j} = [sub_title_str(j);v'];
        end
    end
    sta_re{index_sel} = [sub_re{:}];
end
sta_re = [sta_re{:}]';
if write_sel
    obj_wd.CloseWord()
    xlswrite(xls_fn,sta_re);
end