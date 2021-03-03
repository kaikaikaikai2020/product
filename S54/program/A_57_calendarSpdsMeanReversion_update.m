clear;

key_str = 'A-57';
addpath(genpath(fullfile(pwd,'jplv7')))
% load('inputDataDaily_VX_20120507', 'tday', 'contracts', 'cl');
%load('inputDataDaily_CL_20120813', 'tday', 'contracts', 'cl');
ticker_pool = {'CL','TU','HG','PA'};
for ticker_ID = 1:length(ticker_pool)
    %ticker = 'TU';
    ticker = ticker_pool{ticker_ID};
    sql_str = ['select year(tradeDate)*10000+month(tradeDate)*100+day(tradeDate) as t, ',...
        'ticker,settlePrice from S50.us_futures where mark4 = "%s" and tradeDate>="2006-01-01" order by tradeDate'];
    X = fetchmysql(sprintf(sql_str,ticker),3);
    % T=[1:length(spot)]';
    % isBadData=~isfinite(spot);
    % spot(isBadData)=[];
    % T(isBadData)=[];
    % res=ols(log(spot), [T ones(size(T, 1), 1)]);
    % 
    % fprintf(1, 'Average annualized spot return=%f\n', 252*smartmean(res.beta(1)));
    %f1 = @(x) [x(end-3:end),x(1:end-4)];
    a=table2cell(X(:,2));
    b=cellfun(@(x)  [x(end-3:end),x(1:end-4)],a,'UniformOutput',false);
    X.ticker=b;
    %X.ticker = varfun(@(x) [x(end-3:end),x(1:end-4)],X.ticker);
    X = unstack(X,'settlePrice','ticker');

    cl = table2array(X(:,2:end));
    tday = X.t;
    contracts = X.Properties.VariableNames(2:end);
    % Fitting gamma to forward curve
    gamma=NaN(size(tday));
    for t=1:length(tday)

        FT=cl(t, :)';
        idx=find(isfinite(FT));
        idxDiff=fwdshift(1, idx)-idx; % ensure consecutive months futures
        if (length(idx) >= 5 && all(idxDiff(1:4)==1))
            FT=FT(idx(1:5)); % only uses the nearest 5 contracts
            T=[1:length(FT)]';
    %         scatter(T, log(FT));
            res=ols(log(FT), [T ones(size(T, 1), 1)]);
            gamma(t)=-12*res.beta(1);
        end
    end
    gamma=fillMissingData(gamma);

    % plot(gamma);


    %print -r300 -djpeg fig5_4
    % hold on;

    % fprintf(1, 'Average annualized roll return=%f\n', smartmean(gamma));
    isGoodData=find(isfinite(gamma));
    results=adf(gamma(isGoodData), 0, 1);
    prt(results);

    gammalag=lag(gamma(isGoodData), 1);  
    deltaGamma=gamma(isGoodData)-gammalag;
    deltaGamma(1)=[]; 
    gammalag(1)=[];
    regress_results=ols(deltaGamma, [gammalag ones(size(gammalag))]);
    halflife=-log(2)/regress_results.beta(1);

    fprintf(1, 'halflife=%f days\n', halflife);
    % halflife=36.394034 days


    lookback=round(halflife);
    ma=movingAvg(gamma, lookback);
    mstd=movingStd(gamma, lookback);
    zScore=(gamma-ma)./mstd;

    % linear mean reversion strategy
    isExpireDate=false(size(cl));
    positions=zeros(size(cl));

    isExpireDate=isfinite(cl) & ~isfinite(fwdshift(1, cl));

    holddays=3*21;
    numDaysStart=holddays+10;
    numDaysEnd=10;
    spreadMonth=12; % No. months between far and near contracts.
    for c=1:length(contracts)-spreadMonth
        expireIdx=find(isExpireDate(:, c));
        expireIdx=expireIdx(end); % There may be some missing data earlier on
        if (c==1)
            startIdx=max(1, expireIdx-numDaysStart);
            endIdx=expireIdx-numDaysEnd;
        else % ensure next front month contract doesn't start until current one ends
            myStartIdx=endIdx+1;
            myEndIdx=expireIdx-numDaysEnd;
            if (myEndIdx-myStartIdx >= holddays)
                startIdx=myStartIdx;
                endIdx=myEndIdx;
            else
                startIdx=NaN;
            end
        end

        if (~isempty(expireIdx) && endIdx > startIdx)
            positions(startIdx:endIdx, c)=-1; % Presume we long spread (long back contract, short front contract)
            positions(startIdx:endIdx, c+spreadMonth)=1;
        end
    end



    positions(isnan(zScore), :)=0;
    positions(zScore > 0, :)=-positions(zScore > 0, :);
    % positions(zScore > 1, :)=-positions(zScore > 1, :);

    ret=smartsum(lag(positions).*(cl-lag(cl, 1))./lag(cl, 1), 2)/2;
    ret(isnan(ret))=0;

    idx=find(tday==20080102);
    % idx=1;

    cumret=cumprod(1+ret(idx:end))-1;
     if ~isempty(tday)
        y_re = cumret; 
        tref = cellfun(@num2str,num2cell(tday(idx:end)),'UniformOutput',false);
        tref = datenum(tref,'yyyymmdd');
        tref = cellstr(datestr(tref,'yyyy-mm-dd'));
        %setfigure
        h = figure_S53(y_re,tref,[]);
        title(sprintf('%s-%s',key_str,ticker))
    else
        figure;
        plot(cumret); % Cumulative compounded return
    end
    %plot(cumret); % Cumulative compounded return

    fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+ret(idx:end)).^(252/length(ret(idx:end)))-1, sqrt(252)*mean(ret(idx:end))/std(ret(idx:end)));

    [maxDD maxDDD]=calculateMaxDD(cumret);
    fprintf(1, 'maxDD=%f maxDDD=%i\n', maxDD, maxDDD);
    % APR=0.083406 Sharpe=1.288661
    % maxDD=-0.053222 maxDDD=206
end
