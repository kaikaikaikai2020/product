%{
�뽫������������߶��Ե�Ȩ�صķ�ʽ������������ҿ�һ�²��Ե����ߣ�sharpe���껯�����ÿ�������ͼ��
cadchf,cadjpy,chfjpy,eurcad,eurchf,eurgbp,eurjpy,eurusd,gbpcad,gbpchf,gbpjpy,usdcad,usdchf, audchf,audjpy,euraud,eurnzd

�������Ǹ��ʼ�Ȩ����
%}
clear

symbols_pair = 'cadchf,cadjpy,chfjpy,eurcad,eurchf,eurgbp,eurjpy,eurusd,gbpcad,gbpchf,gbpjpy,usdcad,usdchf,audchf,audjpy,euraud,eurnzd';
symbols_pair = strsplit(symbols_pair,',');
symbols_pair = cellfun(@deblankl,symbols_pair,'UniformOutput',false);

x1 = load('foreign_exchange1.mat');
x2 = load('foreign_exchange2.mat');

symbols = [x1.symbols';x2.symbols];
symbols = cellfun(@deblankl,symbols,'UniformOutput',false);
data = [x1.re;x2.re];

T_pair = length(symbols_pair);
x = cell(T_pair,1);
for i = 1:T_pair
    x{i} = data{strcmp(symbols,symbols_pair(i))}(:,[1:2,end])';
end
t = [x{:}]';
t = unique(t(:,1));

r = zeros(length(t),T_pair);
w = r;
for i = 1:T_pair
    temp = x{i}';
    [~,ia,ib] = intersect(t,temp(:,1));
    r(ia,i) = temp(ib,2);
    w(ia,i) = temp(ib,3);
end

w = bsxfun(@rdivide,w,sum(w,2));
r1 = r.*w;
r_c = cumprod(1+sum(r1,2));
title_str = '���';
t_str = cellstr(datestr(t,'yyyymmdd'));
T = length(t_str);
h1=figure;
subplot(2,1,1)
plot(r_c,'-','LineWidth',2);
set(gca,'xlim',[0,T]);
set(gca,'XTick',floor(linspace(1,T,15)));
set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
set(gca,'XTickLabelRotation',90)    
setpixelposition(h1,[223,365,1345,600]);
box off
title(sprintf('%s-�������',title_str))
%ÿ������
t_y = year(t);
t_y_u = unique(t_y);
r_year = zeros(size(t_y_u));
for j = 1:length(t_y_u)
sub_r = r1(eq(t_y,t_y_u(j)),:);
temp = cumprod(1+sum(sub_r,2));
r_year(j) = temp(end)-1;
end
subplot(2,1,2)
bar(t_y_u,r_year)
box off
title(sprintf('%s-ÿ������ͳ��',title_str))

%�������ߵĲ���
r_1 = [r,sum(r1,2)];
r_str = [symbols_pair,'���'];
sub_re = cell(T_pair,1);
for j = 1:T_pair+1
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

