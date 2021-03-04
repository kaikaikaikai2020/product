%{
请将下面的货币对以等权重的方式结合起来，让我看一下策略的曲线，sharpe和年化收益和每年的收益图形
cadchf,cadjpy,chfjpy,eurcad,eurchf,eurgbp,eurjpy,eurusd,gbpcad,gbpchf,gbpjpy,usdcad,usdchf, audchf,audjpy,euraud,eurnzd
%}
clear

symbols_pair = 'cadchf,cadjpy,eurchf,eurgbp, audchf,audjpy,euraud,eurnzd';
symbols_pair = strsplit(symbols_pair,',');
symbols_pair = [symbols_pair(1:2:end)',symbols_pair(2:2:end)'];
symbols_pair(:,1) = cellfun(@deblankl,symbols_pair(:,1),'UniformOutput',false);
symbols_pair(:,2) = cellfun(@deblankl,symbols_pair(:,2),'UniformOutput',false);

x1 = load('foreign_exchange1.mat');
x2 = load('foreign_exchange2.mat');

symbols = [x1.symbols';x2.symbols];
symbols = cellfun(@deblankl,symbols,'UniformOutput',false);
data = [x1.re;x2.re];

T_pair = size(symbols_pair,1);
sta_re = cell(T_pair,1);

for i = 1:T_pair
    symbol1 = symbols_pair(i,1);
    x = data{strcmp(symbols,symbol1)};

    symbol2 = symbols_pair(i,2);
    y = data{strcmp(symbols,symbol2)};
    t = unique([x(:,1);y(:,1)]);
    r = zeros(length(t),2);
    [~,ia,ib] = intersect(t,x(:,1));
    r(ia,1) = x(ib,2);
    [~,ia,ib] = intersect(t,y(:,1));
    r(ia,2) = y(ib,2);

    r_c = cumprod(1+mean(r,2));
    title_str = strjoin(symbols_pair(i,:),'-');
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
    r_str = [symbols_pair(i,:),strjoin(symbols_pair(i,:),'-')];
    sub_re = cell(3,1);
    for j = 1:3
        sub_y = cumprod(1+r_1(:,j));
        [v0,v_str0] = curve_static(sub_y,[],false);
        [v,v_str] = ad_trans_sta_info(v0,v_str0);
        if eq(j,1)
            sub_re{j} = [[{''},r_str(j)];[v_str;v]'];
        else
            sub_re{j} = [r_str(j);v'];
        end
    end

    sta_re{i} = [sub_re{:}];
end

sta_re = [sta_re{:}]';

