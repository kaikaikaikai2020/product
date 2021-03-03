clear;
addpath(genpath(fullfile(pwd,'jplv7')))
do_sel = 1;
key_str = 'A-51';
if eq(do_sel,0)
    usdcad=load('inputData_USDCAD_20120426', 'tday', 'hhmm', 'cl');
    audusd=load('inputData_AUDUSD_20120426', 'tday', 'hhmm', 'cl');
    firstDate=20090101;
    idx=find(usdcad.tday>firstDate & usdcad.hhmm==1659);
    usdcad.tday=usdcad.tday(idx);
    usdcad.hhmm=usdcad.hhmm(idx);
    usdcad.cl=usdcad.cl(idx);

    idx=find(audusd.tday>firstDate & audusd.hhmm==1659);
    audusd.tday=audusd.tday(idx);
    audusd.hhmm=audusd.hhmm(idx);
    audusd.cl=audusd.cl(idx);
    tday=audusd.tday;
    usdcad_cl=usdcad.cl;
    audusd_cl=audusd.cl;
else
    sql_str = ['select date(tradeDate),closePrice from polygon_fx_minute.%s ',...
        'where tradeDate>="2009-01-01" and 100*hour(tradeDate)+minute(tradeDate)=359 order by tradeDate'];
    usdcad = fetchmysql(sprintf(sql_str,'USDCAD'),2);
    audusd = fetchmysql(sprintf(sql_str,'AUDUSD'),2);
    [tref,ia,ib] = intersect(usdcad(:,1),audusd(:,1));
    usdcad_cl=cell2mat(usdcad(ia,2));
    audusd_cl = cell2mat(audusd(ib,2));
    tday = cellfun(@(x) str2double(strjoin(strsplit(x,'-'),'')),tref);
end
% Need to invert currency pair so that each unit has same capital in local
% currency
cad=1./usdcad_cl;
aud=audusd_cl;

y=[ aud cad   ];
trainlen=250;
lookback=20;
hedgeRatio=NaN(size(y));
numUnits=NaN(size(y, 1), 1);

for t=trainlen+1:size(y, 1)
    res=johansen(y(t-trainlen:t-1, :), 0, 1);
    hedgeRatio(t, :)=res.evec(:, 1)';

    % yport is the market value of a unit portfolio of AUDUSD and CADUSD expressed in US$.
    yport=sum(y(t-lookback+1:t, :).*repmat(hedgeRatio(t, :), [lookback 1]), 2);
    ma=mean(yport);
    mstd=std(yport);
    zScore=(yport(end)-ma)/mstd;

    % numUnits are number of units of unit portfolio of AUDUSD and CADUSD
    numUnits(t)=-(yport(end)-ma)/mstd;

end
tref = tref(trainlen+1:end);

% positions are market values of AUDUSD and CADUSD in portfolio expressed
% in US$.
positions=repmat(numUnits, [1 size(y, 2)]).*hedgeRatio.*y; 

% daily P&L of portfolio in US$.
pnl=sum(lag(positions, 1).*(y-lag(y, 1))./lag(y, 1), 2); 
ret=pnl./sum(abs(lag(positions, 1)), 2); 
ret(isnan(ret))=0;


%plot(cumprod(1+ret(trainlen+1:end))-1); % Cumulative compounded return
if ~isempty(tref)
    y_re = cumprod(1+ret(trainlen+1:end))-1;
    %setfigure
    h = figure_S53(y_re,tref,[]);
    title(sprintf('%s',key_str))
else
    figure;
    plot(cumprod(1+ret(trainlen+1:end))-1); % Cumulative compounded return
end


fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+ret(trainlen+1:end)).^(252/length(ret(trainlen+1:end)))-1, sqrt(252)*mean(ret(trainlen+1:end))/std(ret(trainlen+1:end)));
% APR=0.112410 Sharpe=1.610890


% Kelly leverage
f=mean(ret(trainlen+1:end))/std(ret(trainlen+1:end))^2;
fprintf(1, 'f=%f\n', f);
% f=23.845328

ret=ret(trainlen+1:end);
%save('../Data/AUDCAD_unequal_ret', 'ret');
