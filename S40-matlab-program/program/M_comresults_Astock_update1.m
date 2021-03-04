%{
A股市场组合
%}

clear
% load re_Astock_update1.mat
% load re_Astock_update1_symbols.mat

load re_Astock_update1_moreandless.mat
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
index_name = {'沪深300','中证500','上证50','上证综指'};
sql_str1 = 'select tradeDate,closeIndex from   yuqerdata.yq_index where symbol = ''%s'' order by tradeDate';
sta_re = cell(size(index_name));
for index_sel = 1:length(index_pool)
    title_str = index_name{index_sel};
    sub_x = fetchmysql(sprintf(sql_str1,index_pool{index_sel}),2);
    sub_symbols = yq_methods.get_index_pool(index_pool{index_sel},sub_x{end,1});

    %限制股票池
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
    legend({'组合曲线有手续费','组合曲线无手续费',index_pool{index_sel}},'NumColumns',3,'Location','best')
    box off
    title(sprintf('%s-组合曲线',title_str))
    %每年收益
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
    title(sprintf('%s-每年收益统计',title_str))
    %统计参数
    sub_re = cell(2,1);
    sub_y = [sub_x_c,r_c];
    sub_title_str = {title_str,sprintf('%s组合',title_str)};
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