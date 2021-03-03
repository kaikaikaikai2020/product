%2007F 2007G 是什么数据？
clear;
close all
addpath(genpath(fullfile(pwd,'jplv7')))
key_str = 'A-63';
%{
load('inputData_ETF', 'syms', 'tday', 'cl');

uso=cl(:, strcmp('USO', syms));
xle=cl(:, strcmp('XLE', syms));
tday_ETF=tday;

load('inputDataDaily_CL_20120502', 'tday', 'contracts', 'cl');
%}
%aa = load('inputDataDaily_CL_20120502','contracts','tday');
temp = get_pool_data({'USO','XLE'});
uso=temp.cl(:,1);
xle = temp.cl(:,2);
tday_ETF = temp.tday;

ticker = 'PA';
sql_str = ['select year(tradeDate)*10000+month(tradeDate)*100+day(tradeDate) as t, ',...
    'ticker,settlePrice from S50.us_futures where mark4 = "%s" and tradeDate>="2000-11-20"  ',...
    'and abs(year(tradeDate)-mark2)<=3 order by tradeDate'];
X = fetchmysql(sprintf(sql_str,ticker),3);
a=table2cell(X(:,2));
b=cellfun(@(x)  [x(end-3:end),x(1:end-4)],a,'UniformOutput',false);
X.ticker=b;
%X.ticker = varfun(@(x) [x(end-3:end),x(1:end-4)],X.ticker);
X = unstack(X,'settlePrice','ticker');
cl = table2array(X(:,2:end));
tday = X.t;
contracts = X.Properties.VariableNames(2:end);


% contracts= cellfun(@(x) [x(2:5),x(end)],contracts,'UniformOutput',false);
% [contracts,ia] = intersect(contracts,aa.contracts);
% cl = cl(:,ia);
% 
% [tday,ia] = intersect(tday,aa.tday);
% cl = cl(ia,:);



ratioMatrix=(fwdshift(1, cl')./cl')'; % back/front

ratio=NaN(size(ratioMatrix, 1), 1);
isExpireDate=false(size(ratio));

isExpireDate=isfinite(cl) & ~isfinite(fwdshift(1, cl));

% Define front month as 40 days to 10 days before expiration
numDaysStart=40;
numDaysEnd=10;

for c=1:length(contracts)-1
    expireIdx=find(isExpireDate(:, c),1,'last');
    if (c==1)
        startIdx=expireIdx-numDaysStart;
        endIdx=expireIdx-numDaysEnd;
    else % ensure next front month contract doesn't start until current one ends
        startIdx=max(endIdx+1, expireIdx-numDaysStart);
        endIdx=expireIdx-numDaysEnd;
    end
        
    if (~isempty(expireIdx)) && startIdx>0
        ratio(startIdx:endIdx)=ratioMatrix(startIdx:endIdx, c);
    end
end

[tday,idxA,idxB]=intersect(tday_ETF, tday);
uso=uso(idxA);
xle=xle(idxA);
ratio=ratio(idxB);

positions=zeros(length(tday), 2);

% Contango, negative roll return, buy spot, short future
contango=find(ratio > 1);
positions(contango, :)=repmat([-1 1], [length(contango) 1]);
% Backwardation, positive roll return, short spot, long future
backwardation=find(ratio < 1);
positions(backwardation, :)=repmat([1 -1], [length(backwardation) 1]);

ret=smartsum(lag(positions, 1).*([uso xle]-lag([uso xle], 1))./lag([uso xle], 1), 2)/2;

ret(isnan(ret))=0;

cumret=cumprod(1+ret)-1;

%plot(cumret);
if ~isempty(tday)
    y_re = cumret; 
    tref = cellfun(@num2str,num2cell(tday),'UniformOutput',false);
    tref = datenum(tref,'yyyymmdd');
    tref = cellstr(datestr(tref,'yyyy-mm-dd'));
    %setfigure
    h = figure_S53(y_re,tref,[]);
    title(sprintf('%s-%s',key_str,ticker))
else
    figure;
    plot(cumret); % Cumulative compounded return
end
    
fprintf(1, 'Avg Ann Ret=%7.4f Sharpe ratio=%4.2f \n',252*smartmean(ret), sqrt(252)*smartmean(ret)/smartstd(ret));
fprintf(1, 'APR=%10.4f\n', prod(1+ret).^(252/length(ret))-1);
[maxDD maxDDD]=calculateMaxDD(cumret);
fprintf(1, 'Max DD =%f Max DDD in days=%i\n\n', maxDD, round(maxDDD));
% Avg Ann Ret= 0.1592 Sharpe ratio=1.05 
% APR=    0.1591
% Max DD =-0.192321 Max DDD in days=487


% isContango=zeros(size(ratio));
% isContango(ratio > 1)=1;
% 
% hold on; plot(isContango, 'r'); hold on;