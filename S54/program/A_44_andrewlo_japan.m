%port_trade.m
%³É·Ö¹É
clear;
key_str = 'A-44japan';
addpath(genpath(fullfile(pwd,'jplv7')))
%'../Data/inputDataOHLCDaily_20120424'
%load('inputDataOHLCDaily_stocks_20120424','stocks');
load('index_japan');


index_pool = index_info;
%index_info = {'I000300','I000905','I000001'};
T = length(index_pool);

re2 = cell(T,2);
fee1 = 1/10000;
fee2 = 1/10000;

for index_id = 1:T
    index_sel = index_pool{index_id};
    stocks = index_ticker{index_id};
    stks = get_pool_data_update(stocks,'2005-01-01','japan');
    %load tempstks.mat
    tday = stks.tday;
    cl=stks.cl;
    op=stks.op;
    tref = stks.tref;

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
    %weights(weights<0)=0;
    weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

    %dailyret=smartsum(weights.*(cl-op)./op, 2)./smartsum(abs(weights), 2);
    dailyret=smartsum(weights.*(cl*(1-fee1)-op*(1+fee2))./(op*(1+fee2)), 2)./smartsum(abs(weights), 2);
    
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
    
    pos = [tref',num2cell(weights)];
    pos = cell2table(pos,'VariableNames',['date',sub_info']);
    
    writetable(pos,sprintf('A_44_andrewlo_csi%s.csv',index_sel))
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