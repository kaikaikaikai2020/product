%bollinger.m
clear;
addpath(genpath(fullfile(pwd,'jplv7')))
% 1 minute data on GLD-USO
%{
tref = [];
load('inputData_ETF', 'tday', 'syms', 'cl');
idxG=find(strcmp('GLD', syms));
idxU=find(strcmp('USO', syms));

x=cl(:, idxG);
y=cl(:, idxU);
%}
%load dataGLDUSD
sql_str = 'select tradeDate,closePrice from polygon.usastock_day where ticker = "%s" order by tradeDate';
x=fetchmysql(sprintf(sql_str,'GLD'),2);
y=fetchmysql(sprintf(sql_str,'USO'),2);
[tref,ia,ib] = intersect(x(:,1),y(:,1));
x = cell2mat(x(ia,2));
y = cell2mat(y(ib,2));


lookback=20; % Lookback set arbitrarily short
hedgeRatio=NaN(size(x, 1), 1);
for t=lookback:size(hedgeRatio, 1)
    %regression_result=ols(y(t-lookback+1:t), [x(t-lookback+1:t) ones(lookback, 1)]);    
    %hedgeRatio(t)=regression_result.beta(1);
    regression_result=regress(y(t-lookback+1:t), [x(t-lookback+1:t) ones(lookback, 1)]);   
    hedgeRatio(t)=regression_result(1);
end

y2=[x y];

yport=sum([-hedgeRatio ones(size(hedgeRatio))].*y2, 2); % The net market value of the portfolio is same as the "spread"
hedgeRatio(1:lookback)=[]; % Removed because hedge ratio is indterminate

if ~isempty(tref)
    tref(1:lookback) = [];
end

yport(1:lookback)=[]; 
y2(1:lookback, :)=[];

% Bollinger band strategy
entryZscore=1;
exitZscore=0;

MA=movingAvg(yport, lookback);
MSTD=movingStd(yport, lookback);
zScore=(yport-MA)./MSTD;

longsEntry=zScore < -entryZscore; % a long position means we should buy EWC
longsExit=zScore > -exitZscore;

shortsEntry=zScore > entryZscore;
shortsExit=zScore < exitZscore;

numUnitsLong=NaN(length(yport), 1);
numUnitsShort=NaN(length(yport), 1);

numUnitsLong(1)=0;
numUnitsLong(longsEntry)=1; 
numUnitsLong(longsExit)=0;
numUnitsLong=fillMissingData(numUnitsLong); % fillMissingData can be downloaded from epchan.com/book2. It simply carry forward an existing position from previous day if today's positio is an indeterminate NaN.

numUnitsShort(1)=0;
numUnitsShort(shortsEntry)=-1; 
numUnitsShort(shortsExit)=0;
numUnitsShort=fillMissingData(numUnitsShort);

numUnits=numUnitsLong+numUnitsShort;
positions=repmat(numUnits, [1 size(y2, 2)]).*[-hedgeRatio ones(size(hedgeRatio))].*y2; % [hedgeRatio -ones(size(hedgeRatio))] is the shares allocation, [hedgeRatio -ones(size(hedgeRatio))].*y2 is the dollar capital allocation, while positions is the dollar capital in each ETF.
pnl=sum(lag(positions, 1).*(y2-lag(y2, 1))./lag(y2, 1), 2); % daily P&L of the strategy
ret=pnl./sum(abs(lag(positions, 1)), 2); % return is P&L divided by gross market value of portfolio
ret(isnan(ret))=0;

if ~isempty(tref)
    y_re = cumprod(1+ret)-1;
    %setfigure
    h = figure_S53(y_re,tref,[]);
else
    figure;
    plot(cumprod(1+ret)-1); % Cumulative compounded return
end
fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+ret).^(252/length(ret))-1, sqrt(252)*mean(ret)/std(ret));
% APR=0.178249 Sharpe=0.964673

% Save this for future use
% save('bollinger', 'hedgeRatio', 'MA', 'MSTD');


