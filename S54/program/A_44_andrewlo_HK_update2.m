%port_trade.m
%成分股
%修正开盘价计算收益
%{
1）得到过去10天的开盘竞价移动平均
2)  得到开盘竞价的成交量（volume）
3）过滤出最大成交金额（开盘竞价移动平均*成交量）的N只股票
4）用这N只股票形成股票池，用Andrewlo part2 算法计算收益
%}

clear;
close all
key_str = 'A-44HK-adj-volupdate';
addpath(genpath(fullfile(pwd,'jplv7')))
%'../Data/inputDataOHLCDaily_20120424'
%load('inputDataOHLCDaily_stocks_20120424','stocks');
load('index_HK');
index_pool = index_info;

data00=load('transdata20201011.mat');
data00=data00.datatrans;

dataV = load('HK_vol');
dataV = dataV.data;
symbol_sel = 30;
%index_info = {'I000300','I000905','I000001'};
T = length(index_pool);
%T = 1;
re2 = cell(T,2);
fee1 = 1.2/1000;
fee2 =1.2/1000;
re3 = re2;

for index_id = 1:T
    figure
    index_sel = index_pool{index_id};
    stocks = index_ticker{index_id};
    
    data0=data00.(index_sel);
    sub_V = dataV.(index_sel);
    stks = get_pool_data_update(stocks,'2005-01-01','HK');
    stocks=stks.stocks;
    %load tempstks.mat
    tday = stks.tday;
    cl=stks.cl;
    op=stks.op;
    tref = stks.tref;
    if size(tref,2)>size(tref,1)
        tref = tref';
    end
    %数据对齐
    
    [tref,ia,ib] = intersect(tref,sub_V.tref);
    tday=tday(ia);
    cl=cl(ia,:);
    op = op(ia,:);
    sub_V.tref = sub_V.tref(ib,:);
    sub_V.V = sub_V.V(ib,:);
    
    [stocks,ia,ib] = intersect(stocks,sub_V.symbol);
    cl=cl(:,ia);
    op = op(:,ia);
    sub_V.V=sub_V.V(:,ib);
    %十天移动均值
    op_ma=movmean(op,[9,0]);
    %成交额
    sub_amount = op_ma.*sub_V.V;
    %调整数据
    op_adj=nan(length(data0.tref),size(cl,2));
    [stocks,ia,ib] = intersect(stocks,cellfun(@(x) sprintf('%0.5d',str2double(x(2:end))),data0.stocks,'UniformOutput',false),'stable');
    op_adj(:,ia) = data0.X(:,ib);

    [tref,ia1,ib1] = intersect(tref,data0.tref);    
    op_adj=op_adj(ib1,:);
    op=op(ia1,:);
    cl=cl(ia1,:);
    tday=tday(ia1);
    sub_amount=sub_amount(ia1,:);
    
    ind = isnan(op_adj) | isnan(op);
    op(ind) = nan;
    cl(ind) = nan;
    op_adj(ind) = nan;
    sub_amount(ind) = nan;
    
    %idxStart=find(tday==20070103);
    %idxEnd=find(tday==20111230);
    % cl is a TxN array of closing prices, where T is the number of trading
    % days, and N is the number of stocks in the S&P 500
    ret=(cl-lag(cl, 1))./lag(cl, 1); % daily returns
    ret2 = (cl*(1-fee1)-lag(cl, 1)*(1+fee2))./(lag(cl, 1)*1+fee2); % daily returns real

    marketRet=smartmean(ret, 2); % equal weighted market index return

    weights=-(ret-repmat(marketRet, [1 size(ret, 2)]));
    %基于amount筛选
    for i = 1:size(weights,1)
        sub_a = sub_amount(i,:);
        sub_a(isnan(sub_a)) = 0;
        [~,ia] = sort(sub_a,'descend');
        weights(i,ia(symbol_sel+1:end)) = 0;
    end
    
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
    subplot(2,1,1);
    figure_S53(y_re,tref,sprintf('%s-%s-part1',key_str,index_sel),0);
    %title()

    re2{index_id,1} = y_re;
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
    %基于amount筛选
    for i = 1:size(weights,1)
        sub_a = sub_amount(i,:);
        sub_a(isnan(sub_a)) = 0;
        [~,ia] = sort(sub_a,'descend');
        weights(i,ia(symbol_sel+1:end)) = 0;
    end
    
    weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);
    
    %dailyret=smartsum(weights.*(cl-op)./op, 2)./smartsum(abs(weights), 2);
    %dailyret=smartsum(weights.*(cl*(1-fee1)-op*(1+fee2))./(op*(1+fee2)), 2)./smartsum(abs(weights), 2);
    dailyret=smartsum(weights.*(cl*(1-fee1)-op_adj*(1+fee2))./(op_adj*(1+fee2)), 2)./smartsum(abs(weights), 2);
    n=size(cl,2);
    info4 = cell(n,1);
    for i = 1:n
        temp = [tref,tref,tref,tref,num2cell([op(:,i),cl(:,i),weights(:,i)])];
        temp(:,1)  = {'p2'};
        temp(:,2)={index_sel};
        temp(:,3) = stocks(i);
        ind= ~isnan(op(:,i));
        temp = temp(ind,:);
        info4{i} = temp';
    end
    info4 = [info4{:}];
    
    
    dailyret(isnan(dailyret))=0;

    %plot(cumprod(1+dailyret)-1); % Cumulative compounded return
    subplot(2,1,2);
    y_re = cumprod(1+dailyret);
    %setfigure
    subplot(2,1,2)
    figure_S53(y_re,tref,sprintf('%s-%s-part2',key_str,index_sel),0);

    re2{index_id,2} = y_re;
    fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+dailyret).^(252/length(dailyret))-1, sqrt(252)*mean(dailyret)/std(dailyret));
    
    sub_info = cellfun(@(x) ['A',x],stks.stocks,'UniformOutput',false);
    
    pos = [tref,num2cell(weights)];
    pos = cell2table(pos,'VariableNames',['date',sub_info']);
    
    %writetable(pos,sprintf('%s%s.csv',key_str,index_sel))
    
    re3{index_id}=[info3,info4];
    
    % APR=0.731553 Sharpe=4.713284

    % annavgret=252*smartmean(dailyret, 1)*100
    % 
    % sharpe=sqrt(252)*smartmean(dailyret, 1)/smartstd(dailyret,1) % Sharpe ratio should be about 0.25
    % 
    % % daily pnl with transaction costs deducted
    % onewaytcost=0.0005; % assume 5 basis points
    % 
    % dailyretMinustcost=dailyret - ...
    %     smartsum(abs(weights./cl-backshift(1, weights)./backshift(1, cl)).*backshift(1, cl), 2).*onewaytcost./smartsum(abs(weights), 2); % transaction costs are only incurred when the weights change
    % 
    % annavgretMinustcost=252*smartmean(dailyretMinustcost, 1)*100
    % 
    % sharpeMinustcost=sqrt(252)*smartmean(dailyretMinustcost, 1)/smartstd(dailyretMinustcost, 1) 
    % 
    % % kelly optimal leverage
    % 
    % f=smartmean(dailyretMinustcost, 1)/smartstd(dailyretMinustcost, 1)^2
end

sta_re = curve_static_batch(re2(:,2),index_info);

re3 = [re3{:}]';
[~,ia] = sort(re3(:,4));
re3 = re3(ia,:);
re3 = cell2table(re3,'VariableNames',{'type','code','date','code2','buy','sel','w'});
re3(isnan(re3.w),:) = [];

writetable(re3,sprintf('%s_detail.csv',key_str))

%writetable(pos,sprintf('A_44_andrewlo_csi%s.csv',index_sel))