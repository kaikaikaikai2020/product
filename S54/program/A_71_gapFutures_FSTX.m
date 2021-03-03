clear;
%getTableFromWeb_mod_adair('https://cn.investing.com/indices/euro-stoxx-50-(fstx)-historical-data', 1)
%https://cn.investing.com/indices/eu-stoxx50-historical-data 我实际用的数据
%FSTX指数数据，可以做一个爬取方法
key_str = 'A-71';
addpath(genpath(fullfile(pwd,'jplv7')))

%update data update_A71data
update_A71data();
entryZscore=0.1;
symbol = 'FSTX';
%{
symbol = 'FSTX';
data=load('inputDataOHLCDaily_20120517', 'syms', 'tday', 'op', 'hi', 'lo', 'cl');
idx=find(strcmp(symbol, data.syms));
op=data.op(:, idx);
hi=data.hi(:, idx);
lo=data.lo(:, idx);
cl=data.cl(:, idx);
%}
%[~,~,data]= xlsread('dataA71.xlsx');
data = fetchmysql('select * from S54.A71data',2);
%[~,~,data]= xlsread('a71-2.csv');
%data = data(2:end,:);
tref_num = datenum(data(:,1));
[tref_num,ia] = sort(tref_num);
tref = cellstr(datestr(tref_num,'yyyy-mm-dd'));

%id=tref_num<datenum(2012,5,17);
%ia = ia(id);

data = data(ia,:);
cl = cell2mat(data(:,2));
hi = cell2mat(data(:,4));
lo = cell2mat(data(:,5));
op = cell2mat(data(:,3));



%{
temp = get_pool_data({'FSTX'});
cl=temp.cl;
hi=temp.hi;
lo=temp.lo;
op=temp.op;
tday=temp.tday;
%}
stdretC2C90d=backshift(1, smartMovingStd(calculateReturns(cl, 1), 90));

%longs= op >= backshift(1, hi).*(1+entryZscore*stdretC2C90d);
%shorts=op <= backshift(1, lo).*(1-entryZscore*stdretC2C90d);
%擅自改变了方向 数据相差10倍左右
longs= op <= backshift(1, hi).*(1+entryZscore*stdretC2C90d);
shorts=op >= backshift(1, lo).*(1-entryZscore*stdretC2C90d);


positions=zeros(size(cl));

positions(longs)=1;
positions(shorts)=-1;

ret=positions.*(op-cl)./op;
ret(isnan(ret))=0;

fprintf(1, '%s APR=%10.4f Sharpe=%4.2f\n', symbol, prod(1+ret).^(252/length(ret))-1, mean(ret)*sqrt(252)/std(ret));
% APR=    0.1327 Sharpe=1.44
cumret=cumprod(1+ret)-1; % compounded ROE

if ~isempty(tref)
    y_re = cumret;
    %setfigure
    h = figure_S53(y_re,tref,[]);
    title(sprintf('%s',key_str))
else
    figure;
    plot(cumret); % Cumulative compounded return
end


[maxDD maxDDD]=calculateMaxDD(cumret);
fprintf(1, 'Max DD =%f Max DDD in days=%i\n\n', maxDD, round(maxDDD));

function update_A71data()
    x = getTableFromWeb_mod_adair('https://cn.investing.com/indices/eu-stoxx50-historical-data', 2);
    [x,var_info] = trans_S54_data(x);

    tN = 'S54.A71data';
    t0 = fetchmysql(sprintf('select max(tradeDate) from %s',tN),2);
    tref_num = datenum(x(:,1));
    ind = tref_num>=datenum(t0);
    if any(ind)
        sub_x = x(ind,:);
        %删除最后一组数据
        exemysql(sprintf('delete from %s where tradeDate>="%s"',tN,t0{1}));
        datainsert_adair(tN,var_info,sub_x);
    end
end
