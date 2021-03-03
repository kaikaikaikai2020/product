%port_trade.m
%成分股
%日内买入价格使用开盘前5分钟均价
clear;
key_str = 'A-44csi';
addpath(genpath(fullfile(pwd,'jplv7')))
%'../Data/inputDataOHLCDaily_20120424'
%load('inputDataOHLCDaily_stocks_20120424','stocks');
index_pool = {'000300','000905','000001'};
index_info = index_pool;
data0=load('adj_csi');
data0.tref = cellfun(@char,data0.tref,'UniformOutput',false);
%index_pool = index_info;
%index_info = {'I000300','I000905','I000001'};
T = length(index_pool);

re2 = cell(T,2);
re3 = re2;

fee1 = 1/1000;
fee2 = 1/1000;

for index_id = 1:T
    index_sel = index_pool{index_id};
    stocks = yq_methods.get_index_pool(index_sel,'2005-01-01');
    stks = get_pool_data_update(stocks,'2010-01-01','csi');
    stocks=stks.stocks;
    %load tempstks.mat
    tday = stks.tday;
    cl=stks.cl;
    op=stks.op;
    tref = stks.tref;
    
    op_adj=nan(length(data0.tref),size(cl,2));
    [stocks,ia,ib] = intersect(stocks,data0.symbol,'stable');
    op_adj(:,ia) = data0.X(:,ib);
    
    [tref,ia1,ib1] = intersect(tref,data0.tref);    
    op_adj=op_adj(ib1,:);
    op=op(ia1,:);
    cl=cl(ia1,:);
    tday=tday(ia1);
    
    
    if size(tref,2)>size(tref,1)
        tref = tref';
    end
    
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
    weights(weights<0)=0;
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
    if ~isempty(tref)
        y_re = cumprod(1+dailyret);
        %setfigure
        h = figure_S53(y_re,tref,[]);
        title(sprintf('%s-%s-part1',key_str,index_sel))
    else
        figure;
        plot(cumprod(1+dailyret)-1); % Cumulative compounded return
    end
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

    marketRet=smartmean(ret, 2); % equal weighted market index return

    weights=-(ret-repmat(marketRet, [1 size(ret, 2)])); % weight of a stock is proportional to the negative distance to the market index.
    weights(weights<0)=0;
    weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

    %dailyret=smartsum(weights.*(cl-op)./op, 2)./smartsum(abs(weights), 2);
    %dailyret=smartsum(weights.*(cl*(1-fee1)-op*(1+fee2))./(op*(1+fee2)), 2)./smartsum(abs(weights), 2);
    dailyret=smartsum(weights.*(cl*(1-fee1)-op_adj*(1+fee2))./(op_adj*(1+fee2)), 2)./smartsum(abs(weights), 2);
    
    dailyret(isnan(dailyret))=0;

    %plot(cumprod(1+dailyret)-1); % Cumulative compounded return
    if ~isempty(tref)
        y_re = cumprod(1+dailyret);
        %setfigure
        h = figure_S53(y_re,tref,[]);
        title(sprintf('%s-%s-part2',key_str,index_sel))
    else
        figure;
        plot(cumprod(1+dailyret)-1); % Cumulative compounded return
    end
    re2{index_id,2} = y_re;
    fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+dailyret).^(252/length(dailyret))-1, sqrt(252)*mean(dailyret)/std(dailyret));
    
    sub_info = cellfun(@(x) ['A',x],stks.stocks,'UniformOutput',false);
    
    pos = [tref,num2cell(weights)];
    pos = cell2table(pos,'VariableNames',['date',sub_info']);
    
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
    
    writetable(pos,sprintf('A_44_andrewlo_csi%s.csv',index_sel))
    
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

%{
sta_re = curve_static_batch(re2(:,2),index_info);

re3 = [re3{:}]';
[~,ia] = sort(re3(:,4));
re3 = re3(ia,:);
re3 = cell2table(re3,'VariableNames',{'type','code','date','code2','buy','sel','w'});
re3(re3.w==0,:) = [];
re3(isnan(re3.w),:) = [];
writetable(re3,'A_44_andrewlo_csi_detail.csv')
%}