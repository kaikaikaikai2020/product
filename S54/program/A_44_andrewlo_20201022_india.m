%port_trade.m
%adj-price计算
%补充印度股票结果
clear;
close all

t_cut = {'NIFTY1501';'NIFTY1511';'NIFTY1521';'NIFTY1531';'NIFTY1541';'NIFTY1551';...
    'NIFTY1601';'NIFTY1611';'NIFTY1621';'NIFTY1631';'NIFTY1641';'NIFTY1651';...
    'NIFTY1701';'NIFTY1711';'NIFTY1721';'NIFTY1731';'NIFTY1741';'NIFTY1751';'NIFTY1801';'NIFTY1811';'NIFTY1901'};
t_cut = t_cut(1:end-1);
re_f = cell(size(t_cut));
for t_ind = 1:length(t_cut)


key_str = 'A-andrewlo补充201012india';
addpath(genpath(fullfile(pwd,'jplv7')))
%'../Data/inputDataOHLCDaily_20120424'
%load('inputDataOHLCDaily_stocks_20120424','stocks');
load data20201012.mat
%使用我们计算的数据替换下：
temp = load('transdata20201022');
temp = temp.datatrans;
sub_sel = t_cut{t_ind};

temp = temp.(sub_sel);
data.NIFTY.stocks=temp.stocks;
data.NIFTY.tref = temp.tref;
data.NIFTY.PX_LAST = temp.X;
index_pool = data.index;

%载入数据
%data00=load('transdata20201011.mat'); %5分钟修正
data00=load('transdata20201016.mat'); %10分钟修正

data00 = data00.datatrans;
data00.KRX100=data00.KRX;
data00.MSCITW=data00.MSCI;
%数据对齐
data00.index=[data00.index;'KRX100';'MSCITW'];
temp_index_pool = cell(size(index_pool));
for i = 1:length(temp_index_pool)
    temp = index_pool{i};
    temp = strsplit(temp,'_');
    temp_index_pool{i} = temp{end};
end
[temp_index_pool,ia] = intersect(temp_index_pool,data00.index);
index_pool = index_pool(ia);
index_info = index_pool;
T = length(index_pool);
re2 = cell(T,2);
re3 = re2;

for index_id =  1:T
    %figure
    index_sel = index_pool{index_id};
    
    temp= index_info{index_id};
    t_str = strrep(temp,'_','-');
    
    temp2 = temp_index_pool{index_id};
    
    sub_data = data.(index_pool{index_id});
    %stocks = sub_data.stocks;
    stocks = sub_data.stocks;
    var_info = fieldnames(sub_data);
    cl=sub_data.(var_info{1});
    op=sub_data.(var_info{2});    
    tref = sub_data.(var_info{7});
    
    data0=data00.(temp2);
    op_adj=nan(length(data0.tref),size(cl,2));
    temp_stocks = cell(size(stocks));
    for i = 1:length(stocks)
        temp =strsplit(stocks{i},' ');
        temp_stocks(i) = temp(1);
    end
    if any(index_id ==[1,2,5])
        temp_stocks1 = data0.stocks;
    elseif index_id ==3
        temp_stocks1 = cellfun(@(x) x(2:end),data0.stocks,'UniformOutput',false);
    else
        temp_stocks1 = cellfun(@(x) x(3:end),data0.stocks,'UniformOutput',false);
    end
    [~,ia,ib] = intersect(temp_stocks,temp_stocks1,'stable');
    op_adj(:,ia) = data0.X(:,ib);

    [tref,ia1,ib1] = intersect(tref,data0.tref);    
    op_adj=op_adj(ib1,:);
    op=op(ia1,:);
    cl=cl(ia1,:);
    
    ind = ~all(isnan(cl'))' & ~all(isnan(op_adj'))';
    cl=cl(ind,:);
    op=op(ind,:);
    op_adj=op_adj(ind,:);
    tref = tref(ind);
    
    tday= cellfun(@(x) str2double(strrep(x,'-','')),tref);
    fee1 = sub_data.fee1;
    fee2 = sub_data.fee2;
    

    %idxStart=find(tday==20070103);
    %idxEnd=find(tday==20111230);
    ind = tday>=20100101;
    %ind = tday>=20070103 & tday<=20111230;
    tday=tday(ind);
    cl=cl(ind, :);
    op=op(ind, :);
    tref=tref(ind);
    % cl is a TxN array of closing prices, where T is the number of trading
    % days, and N is the number of stocks in the S&P 500
    ret=(cl-lag(cl, 1))./lag(cl, 1); % daily returns
    ret2 = (cl*(1-fee1)-lag(cl, 1)*(1+fee2))./(lag(cl, 1)*1+fee2); % daily returns real

    marketRet=smartmean(ret, 2); % equal weighted market index return

    weights=-(ret-repmat(marketRet, [1 size(ret, 2)]));
    %weights(weights<0)=0;
    weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

    %dailyret=smartsum(backshift(1, weights).*ret, 2); % Capital is always one
    dailyret=smartsum(backshift(1, weights).*ret2, 2); % Capital is always one
    dailyret(isnan(dailyret))=0;
    
    n=size(cl,2);
    info3 = cell(n,1);
    for i = 1:n
        temp = [tref,tref,tref,tref,num2cell([lag(cl(:,i), 1),cl(:,i),backshift(1, weights(:,i))])];
        temp(:,1)  = {'p1'};
        temp(:,2)={index_sel};
        temp(:,3) = stocks(i);
        ind= ~isnan(lag(cl(:,i),1)) & ~eq(lag(cl(:,i),1),0);
        temp = temp(ind,:);
        info3{i} = temp';
    end
    info3 = [info3{:}];
    
    %plot(cumprod(1+dailyret)-1); % Cumulative compounded return
    y_re = cumprod(1+dailyret);
    %setfigure
    %temp=sprintf('%s-%s-part1',key_str,t_str);
    %subplot(2,1,1);
    %h = figure_S53(y_re,tref,temp,0);        

    re2{index_id,1} = y_re;
    fprintf(1,'%s 隔日交易\n',t_str);
    fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+dailyret).^(252/length(dailyret))-1, sqrt(252)*mean(dailyret)/std(dailyret));
    % APR=13.7%, Sharpe=1.3

    % daily pnl with transaction costs deducted
    % onewaytcost=0.0005; % assume 5 basis points
    % 
    % dailyretMinustcost=dailyret - ...
    %     smartsum(abs(weights./cl-backshift(1, weights)./backshift(1, cl)).*backshift(1, cl), 2).*onewaytcost./smartsum(abs(weights), 2); % transaction costs are only incurred when the weights change
    % 
    % annavgretMinustcost=252*smartmean(dailyretMinustcost, 1)*100
    % 
    % sharpeMinustcost=sqrt(252)*smartmean(dailyretMinustcost, 1)/smartstd(dailyretMinustcost, 1) 
    % 
    % % switch to use open prices
    % 
    ret=(op-backshift(1, cl))./backshift(1, cl); % daily returns
    %ret2=(op*(1-fee1)-backshift(1, cl)*(1+fee2))./(backshift(1, cl)*(1+fee2); % daily returns

    marketRet=smartmean(ret, 2); % equal weighted market index return

    weights=-(ret-repmat(marketRet, [1 size(ret, 2)])); % weight of a stock is proportional to the negative distance to the market index.
    %weights(weights<0)=0;
    weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

    %dailyret=smartsum(weights.*(cl-op)./op, 2)./smartsum(abs(weights), 2);
    dailyret=smartsum(weights.*(cl*(1-fee1)-op_adj*(1+fee2))./(op_adj*(1+fee2)), 2)./smartsum(abs(weights), 2);
    dailyret0=smartsum(weights.*(cl*(1-fee1)-op*(1+fee2))./(op*(1+fee2)), 2)./smartsum(abs(weights), 2);
    
    n=size(cl,2);
    info4 = cell(n,1);
    for i = 1:n
        temp = [tref,tref,tref,tref,num2cell([op_adj(:,i),cl(:,i),weights(:,i)])];
        temp(:,1)  = {'p2'};
        temp(:,2)={index_sel};
        temp(:,3) = stocks(i);
        ind= ~isnan(op(:,i));
        temp = temp(ind,:);
        info4{i} = temp';
    end
    info4 = [info4{:}];
    
    
    dailyret(isnan(dailyret))=0;
    dailyret0(isnan(dailyret0)) = 0;
    
    %plot(cumprod(1+dailyret)-1); % Cumulative compounded return
    y_re0 = cumprod(1+dailyret0);
    y_re = cumprod(1+dailyret);
    %setfigure
    temp=sprintf('%s-%s-part2',key_str,t_str);
    %subplot(2,1,2)
    %h = figure_S53([y_re],tref,temp,0);
    re_f{t_ind} = y_re;
    %legend({'修正后'})
    re2{index_id,2} = y_re;
    fprintf(1,'%s 日内交易\n',t_str);
    fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+dailyret).^(252/length(dailyret))-1, sqrt(252)*mean(dailyret)/std(dailyret));
end

end


y = re_f;
y = [y{:}];
h = figure_S53(y,tref,'NIFTY时间测试',1);
legend(t_cut,'location','bestoutside')
