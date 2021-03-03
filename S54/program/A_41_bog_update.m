clear;
%大约1分钟
%难点 S&P 500成分股代号获取
addpath(genpath(fullfile(pwd,'jplv7')))
topN=10; % Max number of positions
entryZscore=1;
lookback=20; % for MA
%load('../Data/inputDataOHLCDaily_20120424', 'stocks', 'tday', 'op', 'hi', 'lo', 'cl');
%load('inputDataOHLCDaily_stocks_20120424', 'stocks', 'tday', 'op', 'hi', 'lo', 'cl');
key_str = 'A_41';
sql_str = 'select * from S54.A41data order by tradeDate';
X = fetchmysql(sql_str,2);
%tref = cellfun(@(x) x(:,1)',X,'UniformOutput',false);
tref = unique(X(:,2));
tday = cellfun(@(x) str2double(strjoin(strsplit(x,'-'))),tref);
stocks = unique(X(:,1));
T = length(stocks);
x=cell(T,1);

parfor i = 1:T
    ind = strcmp(X(:,1),stocks(i));
    temp = X(ind,2:end);
    if ~isempty(temp)
        [~,ia,ib] = intersect(tref,temp(:,1));
        sub_x = nan(length(tref),size(temp,2)-1);
        sub_x(ia,:) = cell2mat(temp(ib,2:end));
    else
        keyboard
    end
    x{i} = sub_x;
end

temp = cell(4,1);
for i = 1:4
    sub_temp = cellfun(@(x) x(:,i),x,'UniformOutput',false);
    temp{i} = [sub_temp{:}];
end

op=temp{1};
hi=temp{2};
lo=temp{3};
cl=temp{4};

stdretC2C90d=backshift(1, smartMovingStd(calculateReturns(cl, 1), 90));
buyPrice=backshift(1, lo).*(1-entryZscore*stdretC2C90d);

retGap=(op-backshift(1, lo))./backshift(1, lo);

pnl=zeros(length(tday), 1);

positionTable=zeros(size(cl));

ma=backshift(1, smartMovingAvg(cl, lookback));

for t=2:size(cl, 1)
    hasData=find(isfinite(retGap(t, :)) & op(t, :) < buyPrice(t, :) & op(t, :) > ma(t, :));
    
    [foo idxSort]=sort(retGap(t, hasData), 'ascend');
    positionTable(t, hasData(idxSort(1:min(topN, length(idxSort)))))=1;
end

retO2C=(cl-op)./op;


pnl=smartsum(positionTable.*(retO2C), 2);
ret=pnl/topN; 
ret(isnan(ret))=0;

fprintf(1, '%i - %i\n', tday(1), tday(end));
fprintf(1, 'APR=%10.4f\n', prod(1+ret).^(252/length(ret))-1);

fprintf(1, 'Sharpe=%4.2f\n', mean(ret)*sqrt(252)/std(ret));
% APR=8.7%, Sharpe=1.5

%cumret=cumprod(1+ret)-1; % compounded ROE

%plot(cumret);
if ~isempty(tref)
    y_re = cumprod(1+ret)-1;
    %setfigure
    h = figure_S53(y_re,tref,[]);
    title('A-41')
else
    figure;
    plot(cumprod(1+ret)-1); % Cumulative compounded return
end
