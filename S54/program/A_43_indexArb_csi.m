%成分股获取
%上证50和沪深300
clear;
key_str = 'A-43';
addpath(genpath(fullfile(pwd,'jplv7')))

index_pool = {'000016','000300'};
etf_pool = {'510050','510300'};
t0 = '2007-01-01';

sql_str = ['select tradeDate,openPrice,HighestPrice,lowestPrice,closePrice ',...
    'from yuqerdata.MktFunddGet where ticker="%s" and tradeDate >="%s" order by tradeDate'];
for index_ID = 1%1:2
    etf_sel = etf_pool{index_ID};    
    etf0 = fetchmysql(sprintf(sql_str,etf_sel,t0),2);
    etf = [];
    etf.cl = cell2mat(etf0(:,5));
    etf.tday = cellfun(@(x) str2double(strjoin(strsplit(x,'-'),'')),etf0(:,1));
    index_sel = index_pool{index_ID};
    stocks = yq_methods.get_index_pool(index_sel,t0);
    stks = get_pool_data_update(stocks,t0,'csi');
    
    %}

    % Ensure data have same dates
    [tday,idx1,idx2]=intersect(stks.tday, etf.tday);

    tref = stks.tref(idx1);

    stks.cl=stks.cl(idx1, :);
    etf.cl=etf.cl(idx2, :);

    % Use SPY
    %idxS=find(strcmp('SPY', etf.syms));
    %etf.cl=etf.cl(:, idxS);

    trainDataIdx=find(tday>=20120101 & tday<=20151231);
    testDataIdx=find(tday > 20151231);
    if ~isempty(tref)
        tref = tref(tday>20151231);
    end

    isCoint=false(size(stks.stocks));
    for s=1:length(stks.stocks)
        % Combine the two time series into a matrix y2 for input into Johansen test
        y2=[stks.cl(trainDataIdx, s), etf.cl(trainDataIdx)];
        badData=any(isnan(y2), 2);
        y2(badData, :)=[]; % remove any missing data

        if (size(y2, 1) > 250)
            results=johansen(y2, 0, 1); % johansen test with non-zero offset but zero drift, and with the lag k=1.
            if (results.lr1(1) > results.cvt(1, 1))
                isCoint(s)=true;
            end
        end    
    end

    length(find(isCoint))
    % 98: there are 98 stocks that are cointegrating with SPY

    % Form a long-only portfolio with all stocks that cointegrate with SPY, with equal
    % capital allocation
    yN=stks.cl(trainDataIdx, isCoint);
    yN=fillmissing(yN,'previous');
    test_id = sum(yN);
    yN(:,isnan(test_id)) = [];
    isCoint(isnan(test_id)) = [];
    logMktVal_long=sum(log(yN), 2); % The net market value of the long-only portfolio is same as the "spread"

    % Confirm that the portfolio cointegrates with SPY
    ytest=[logMktVal_long, log(etf.cl(trainDataIdx))]; 
    results=johansen(ytest, 0, 1); % johansen test with non-zero offset but zero drift, and with the lag k=1.
    prt(results);

    % Output:
    %  Johansen MLE estimates 
    % NULL:                  Trace Statistic         Crit 90%         Crit 95%         Crit 99% 
    % r <= 0   variable   1           15.869           13.429           15.494           19.935 
    % r <= 1   variable   2            6.197            2.705            3.841            6.635 
    % 
    % NULL:                  Eigen Statistic         Crit 90%         Crit 95%         Crit 99% 
    % r <= 0   variable   1            9.671           12.297           14.264           18.520 
    % r <= 1   variable   2            6.197            2.705            3.841            6.635 

    results.evec
    % 
    % ans =
    % 
    %     1.0939   -0.2799
    %  -105.5600   56.0933


    % Apply linear mean-reversion model on test set
    yNplus=[stks.cl(testDataIdx, isCoint), etf.cl(testDataIdx)]; % Array of stock and ETF prices
    weights=[repmat(results.evec(1, 1), size(stks.cl(testDataIdx, isCoint))), ...
           repmat(results.evec(2, 1), size(etf.cl(testDataIdx)))]; % Array of log market value of stocks and ETF's

    logMktVal=smartsum(weights.*log(yNplus), 2); % Log market value of long-short portfolio

    lookback=5;
    numUnits=-(logMktVal-movingAvg(logMktVal, lookback))./movingStd(logMktVal, lookback); % capital invested in portfolio in dollars.  movingAvg and movingStd are functions from epchan.com/book2
    positions=repmat(numUnits, [1 size(weights, 2)]).*weights; % positions is the dollar capital in each stock or ETF.
    pnl=smartsum(lag(positions, 1).*(log(yNplus)-lag(log(yNplus), 1)), 2); % daily P&L of the strategy
    ret=pnl./smartsum(abs(lag(positions, 1)), 2); % return is P&L divided by gross market value of portfolio
    ret(isnan(ret))=0;

    %figure;
    %plot(cumprod(1+ret)-1); % Cumulative compounded return
    if ~isempty(tref)
        y_re = cumprod(1+ret)-1;
        %setfigure
        h = figure_S53(y_re,tref,[]);
        title(sprintf('%s-%s-%s',key_str,index_sel,etf_sel))
    else
        figure;
        plot(cumprod(1+ret)-1); % Cumulative compounded return
    end

    fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+ret).^(252/length(ret))-1, sqrt(252)*mean(ret)/std(ret));
    % APR=0.044930 Sharpe=1.319397
end