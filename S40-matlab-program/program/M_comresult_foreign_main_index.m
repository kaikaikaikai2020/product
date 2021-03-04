%{
请将下面的所有曲线对以等权重的方式结合起来，让我看一下策略的曲线，sharpe和年化收益和每年的收益图形
cadchf,cadjpy,chfjpy,eurcad,eurchf,eurgbp,eurjpy,eurusd,gbpcad,gbpchf,gbpjpy,usdcad,usdchf, audchf,audjpy,euraud,eurnzd
%}
clear

load foreign_main_index_update2

del_ind = cellfun(@isempty,re);
re(del_ind) = [];
symbols(del_ind) = [];
symbols = symbols';

symbols_pair = symbols;
T_pair = length(symbols_pair);

x = cell(size(symbols_pair));
for i = 1:length(x)
    if ~isempty(re{i})
        x{i} = re{i}(:,1:2)';
    end
end



t = [x{:}]';
t = unique(t(:,1));

r = zeros(length(t),T_pair);
for i = 1:T_pair
    temp = x{i}';
    [~,ia,ib] = intersect(t,temp(:,1));
    r(ia,i) = temp(ib,2);
end
r_c = cumprod(1+mean(r,2));
title_str = '组合';
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
title(sprintf('%s-组合曲线',title_str))
%每年收益
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
title(sprintf('%s-每年收益统计',title_str))

%三条曲线的参数
r_1 = [r,mean(r,2)];
r_str = [symbols_pair,'组合'];
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

