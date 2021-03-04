%{
���޼��������
audcad,audchf,audjpy,audnzd,audusd,nzdusd,nzdjpy,euraud,eurnzd,usdjpy
2:00am to 8:00am �Ǽ������ʱ���
8.00am to 14:00 �ǽ���ʱ���

���޼����������ʣ���
audcad,audchf,audjpy,audnzd,audusd,nzdusd,nzdjpy,euraud,eurnzd,usdjpy
8:00am to 14:00 �Ǽ������ʱ���
14.00am to 20:00 �ǽ���ʱ���
%}

clear

write_sel = true;
key_str = 'S40SMT������Ҫ���';

dN = 'foreign_index_min_V2';
symbols = strsplit('cadchf,cadjpy,chfjpy,eurcad,eurchf,eurgbp,eurjpy,eurusd,gbpcad,gbpchf,gbpjpy,usdcad,usdchf,audchf,audjpy,euraud,eurnzd',',');
symbols1 = strsplit('audcad,audchf,audjpy,audnzd,audusd,nzdusd,nzdjpy,euraud,eurnzd,usdjpy',',');
T_symbols = length(symbols);
re = cell(T_symbols,1);
X = re;
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
    P.out_sel = false;
    
    title_str = symbols{i_sym};  
    if any(strcmp(symbols1,title_str))
        t_num1 = 200;
        t_num2 = 1400;
    else
        t_num1 = 800;
        t_num2 = 2000;
    end
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
    id_sel = temp>t_num1&temp<=t_num2;
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
    [tradeYield,result1,tradeDetail,yearDetail,h,assure_ratio] = SMTTradingModelTool2(closeprice,openprice,t,closeprice,openprice,t,P);
    re{i_sym} = [tradeYield,assure_ratio];
    X{i_sym} = tradeYield(:,1:2)';
    sprintf('%s %d-%d',key_str,i_sym,T_symbols)
    
end

t = [X{:}]';
t = unique(t(:,1));

r = zeros(length(t),T_symbols);
for i = 1:T_symbols
    temp = X{i}';
    [~,ia,ib] = intersect(t,temp(:,1));
    r(ia,i) = temp(ib,2);
end
r_c = cumprod(1+mean(r,2));
title_str = '���';
t_str = cellstr(datestr(t,'yyyymmdd'));
T = length(t_str);
h=figure;
subplot(2,1,1)
plot(r_c,'-','LineWidth',2);
set(gca,'xlim',[0,T]);
set(gca,'XTick',floor(linspace(1,T,15)));
set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
set(gca,'XTickLabelRotation',90)    
setpixelposition(h,[223,365,1345,600]);
box off
title(sprintf('%s-�������',title_str))
%ÿ������
t_y = year(t);
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

%�������ߵĲ���
r_1 = [r,mean(r,2)];
r_str = [symbols,'���'];
sub_re = cell(T_symbols,1);
for j = 1:T_symbols+1
    sub_y = cumprod(1+r_1(:,j));
    [v0,v_str0] = curve_static(sub_y,[],false);
    [v,v_str] = ad_trans_sta_info(v0,v_str0);
    if eq(j,1)
        sub_re{j} = [[{''},r_str(j)];[v_str;v]'];
    else
        sub_re{j} = [r_str(j);v'];
    end
end
sub_re = [sub_re{:}]';

if write_sel
    pn_write = fullfile(pwd,'������');
    if ~exist(pn_write,'dir')
        mkdir(pn_write)
    end
    obj_wd = wordcom(fullfile(pn_write,sprintf('%s.doc',key_str)));
    xls_fn = fullfile(pn_write,sprintf('%s.xlsx',key_str));
    
    obj_wd.pasteFigure(h,key_str);
    obj_wd.CloseWord();
    xlswrite(xls_fn,sub_re);
end