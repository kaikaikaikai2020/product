%port_trade.m
%³É·Ö¹É
clear;
key_str = 'A-44';
addpath(genpath(fullfile(pwd,'jplv7')))
%'../Data/inputDataOHLCDaily_20120424'
%load('inputDataOHLCDaily_stocks_20120424','stocks');
%stks=load('inputDataOHLCDaily_stocks_20120424');
%stks.tref = [];
stks = get_A41data();
sub_info = cellfun(@(x) replace(x,'.','_'),stks.stocks,'UniformOutput',false);
%stks = get_pool_data_update(stocks);
tday = stks.tday;
cl=stks.cl;
op=stks.op;
tref = stks.tref;

%idxStart=find(tday==20070103);
%idxEnd=find(tday==20111230);
ind = tday>=20070103;
%ind = tday>=20070103 & tday<=20111230;
tday=tday(ind);
cl=cl(ind, :);
op=op(ind, :);
if ~isempty(tref)
tref=tref(ind);
end
% cl is a TxN array of closing prices, where T is the number of trading
% days, and N is the number of stocks in the S&P 500
ret=(cl-lag(cl, 1))./lag(cl, 1); % daily returns

marketRet=smartmean(ret, 2); % equal weighted market index return

weights=-(ret-repmat(marketRet, [1 size(ret, 2)]));
weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

dailyret=smartsum(backshift(1, weights).*ret, 2); % Capital is always one

dailyret(isnan(dailyret))=0;

%plot(cumprod(1+dailyret)-1); % Cumulative compounded return
if ~isempty(tref)
    y_re = cumprod(1+dailyret)-1;
    %setfigure
    h = figure_S53(y_re,tref,[]);
    title(sprintf('%s-part1',key_str))
else
    figure;
    plot(cumprod(1+dailyret)-1); % Cumulative compounded return
end
pos = [tref,num2cell(weights)];
pos = cell2table(pos,'VariableNames',['date',sub_info']);
writetable(pos,'A_44_andrewlo_2007_2012_closeprice.csv')

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
weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

dailyret=smartsum(weights.*(cl-op)./op, 2)./smartsum(abs(weights), 2);
dailyret(isnan(dailyret))=0;

%plot(cumprod(1+dailyret)-1); % Cumulative compounded return
if ~isempty(tref)
    y_re = cumprod(1+dailyret)-1;
    %setfigure
    h = figure_S53(y_re,tref,[]);
    title(sprintf('%s-part2',key_str))
else
    figure;
    plot(cumprod(1+dailyret)-1); % Cumulative compounded return
end
fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+dailyret).^(252/length(dailyret))-1, sqrt(252)*mean(dailyret)/std(dailyret));

pos = [tref,num2cell(weights)];
pos = cell2table(pos,'VariableNames',['date',sub_info']);
writetable(pos,'A_44_andrewlo_2007_2012_openprice.csv')

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
