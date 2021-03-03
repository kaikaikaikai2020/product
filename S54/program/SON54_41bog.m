%国内成分计算
%日本指数（TPX-TOPIX无数据   NKY-N225）、香港指数（HSI、HSCEI）、美国指数（SPY）csi300 500
clear;
%难点 S&P 500成分股代号获取

key_str = 'A41Bog综合计算';
path_save='计算结果';
if ~exist(path_save,'dir')
    mkdir(path_save)
end
addpath(genpath(fullfile(pwd,'jplv7')))
topN=10; % Max number of positions
entryZscore=1;
lookback=20; % for MA

load A41Bog_data.mat
%load('bloombergdata.mat');
index_pool = index;
index1 = cellfun(@(x) ['I_',x],index,'UniformOutput',false);
%index_info = {'I000300','I000905','I000001'};
T = length(index_pool);
%T=2;
re = cell(T,1);
tref0 = [];
re2 =cell(T,1);
re3 = cell(size(re2));
re4 = cell(T,1);
re5 = cell(T,1);
re6 = cell(T,1);
re7 = cell(T,1); % 交易信号中间结果

for index_id = 1:T
    stocks = index_com{index_id};
    dtype=index_info{index_id};
    if strcmpi(dtype,'HK')
        fee1 = 1.2/1000;
        fee2 =1.2/1000;
    elseif strcmpi(dtype,'japan')
        fee1 = 1/10000;
        fee2 = 1/10000;
    elseif strcmpi(dtype,'us')
        fee1 = 1/10000;
        fee2 = 1/10000;
    elseif strcmpi(dtype,'csi')
        fee1 = 1/1000;
        fee2 = 1/1000;
    else
        sprintf('%s 未知品种',key_str)
        fee1 = 1/1000;
        fee2 = 1/1000;
    end
    temp = get_pool_data(stocks,'2012-01-01',dtype);
    stocks = temp.stocks;
    cl=temp.cl;
    hi=temp.hi;
    lo=temp.lo;
    op=temp.op;
    tday=temp.tday;
    tref = temp.tref;

    stdretC2C90d=backshift(1, smartMovingStd(calculateReturns(cl, 1), 90));    
    stdretC2C90d0 = smartMovingStd(calculateReturns(cl, 1), 90);    
    buyPrice=backshift(1, lo).*(1-entryZscore*stdretC2C90d);    
    buyPrice0=lo.*(1-entryZscore*stdretC2C90d0);
    retGap=(op-backshift(1, lo))./backshift(1, lo);
    %pnl=zeros(length(tday), 1);

    positionTable=zeros(size(cl));

    ma=backshift(1, smartMovingAvg(cl, lookback));
    ma0=smartMovingAvg(cl, lookback);

    for t=2:size(cl, 1)
        hasData=find(isfinite(retGap(t, :)) & op(t, :) < buyPrice(t, :) & op(t, :) > ma(t, :));

        [foo idxSort]=sort(retGap(t, hasData), 'ascend');
        positionTable(t, hasData(idxSort(1:min(topN, length(idxSort)))))=1;
    end
    
    temp1 =[stocks,stocks,stocks, num2cell([buyPrice0(end,:)',ma0(end,:)',lo(end,:)'])];  
    temp1(:,2) = index_pool(index_id);
    temp1(:,1) = {datestr(datenum(tref(end))+1,'yyyy-mm-dd')};
    %temp1 = cell2table(temp1,'VariableNames',{'tradeDate','Index','com_ticker','buyPrice','ma','lo'});
    
    re7{index_id} = temp1';
    
    %retO2C=(cl-op)./op;
    retO2C=(cl*(1-fee1)-op*(1+fee2))./(op*(1+fee2));
    pnl=smartsum(positionTable.*(retO2C), 2);
    %fee2 = smartsum(positionTable.*(fee*ones(size(retO2C))),2);
    
    
    ret=pnl/topN; 
    ret(isnan(ret))=0;    
    %fee2 = fee2/topN;
    %ret = ret-fee2;    
    sub_pool = cell(size(positionTable(:,1)));
    sub_pool_num = zeros(size(positionTable(:,1)));
    for i = 1:length(sub_pool)
        temp = positionTable(i,:);
        sub_sym = stocks(temp>0);
        sub_buy = op(i,temp>0);
        sub_close=cl(i,temp>0);
        sub_pool_num(i) = length(sub_sym);
        if ~isempty(sub_sym)
            sub_pool{i} = strjoin(sub_sym,'-');
        else
            sub_pool{i} = '_';
        end
    end
    
    fprintf(1, '%i - %i\n', tday(1), tday(end));
    fprintf(1, 'APR=%10.4f\n', prod(1+ret).^(252/length(ret))-1);
    fprintf(1, 'Sharpe=%4.2f\n', mean(ret)*sqrt(252)/std(ret));
    %cumret=cumprod(1+ret)-1; % compounded ROE
    if size(tref,2)>1
        tref = tref';
    end
    %plot(cumret);
    if ~isempty(tref)
        y_re = cumprod(1+ret)-1;
        %setfigure
        figure;
        subplot(2,1,1)
        figure_S53(y_re,tref,sprintf('A-41-%s',index{index_id}),0);
    else
        y_re = cumprod(1+ret)-1;
        figure;
        plot(y_re); % Cumulative compounded return
    end
    
    pos0=[tref,sub_pool];
    pos0=cell2table(pos0,'VariableNames',{'date',index1{index_id}});
    tref0=unique([tref0;table2cell(pos0(:,1))]);
    re{index_id}=pos0;
    re2{index_id} = y_re;
    %%{
    info2 = cell(size(positionTable(:,1)));
    for i =  1:length(sub_pool)
        temp = positionTable(i,:);
        sub_sym = stocks(temp>0);
        sub_buy = op(i,temp>0);
        sub_close=cl(i,temp>0);
        if ~isempty(sub_sym)
            temp = [sub_sym,num2cell([sub_buy;sub_close]')];
            temp = temp(:,[1,1,1:end]);
            temp(:,2) = tref(i);
            temp(:,1) = index(index_id);
            info2{i} = temp';
        end
    end
    info2 = [info2{:}];
    re3{index_id} = info2;
    %}
    
    temp = get_year_return(tref,ret);
    temp(:,1) = index_pool(index_id);
    
    re4{index_id} = temp';
    
    %%{
    %指数对冲部分，暂时不接入
    db_sel=index_db{index_id};
    sub_ticker = index_0{index_id};
    
    if strcmpi(db_sel,'uqer')
        tN='yuqerdata.yq_index';
        sql_tmp = 'select tradeDate,closeIndex from %s where symbol = "%s" order by tradeDate';
    else
        tN='S50.index_sina';
        sql_tmp = 'select tradeDate,closeIndex from %s where ticker = "%s" order by tradeDate';
    end
    
    sub_index=fetchmysql(sprintf(sql_tmp,tN,sub_ticker),2);
    %disp(size(x));
    %index=data.(index_pool{index_id});
    [tref2,ia,ib] = intersect(tref,sub_index(:,1));
    sub_index = cell2mat(sub_index(ib,2));
    sub_index(2:end) = sub_index(2:end)./sub_index(1:end-1)-1;
    sub_index(1) = 0;
    
    ret=ret(ia,:);
    ret(1) = 0;
    sub_pool_num = sub_pool_num(ia);
    ret_ls=(ret-sub_index.*sub_pool_num)./(1+sub_pool_num);
    y_ls=cumprod(1+ret_ls);    
    %对冲结果
    
    temp = get_year_return(tref2,ret_ls);
    temp(:,1) = {sprintf('%sls',index_pool{index_id})};    
    re5{index_id} = temp';
    re6{index_id} = y_ls;
    subplot(2,1,2)
    figure_S53(y_ls-1,tref2,sprintf('A41bog-%s-LS',index_pool{index_id}),0);
    %}
end

%%{
re3 = [re3{:}]';
re3 = cell2table(re3,'VariableNames',{'code','date','code2','buy','sel'});
%writetable(re3,'A_41_bog_HK_detail.csv')
fn00=fullfile(path_save,sprintf('%s-历史交易明细.csv',key_str));
writetable(re3,fn00)

%}
X = cell(T,1);
var =cell(T,1);
for i = 1:T
    sub_x = re{i};
    temp = sub_x.Properties.VariableNames(2:end);
    temp = cellfun(@(x) sprintf('%s_I%d',x,i),temp,'UniformOutput',false);
    var{i} = temp;
    sub_x = table2cell(sub_x);
    sub_x2 = cell(length(tref0),size(sub_x,2)-1);
    [~,ia,ib] = intersect(tref0,sub_x(:,1),'stable');
    sub_x2(ia,:) = sub_x(ib,2:end);
    X{i} = sub_x2;    
end

X=[tref0,[X{:}]];
var=['date',[var{:}]];
X = cell2table(X,'VariableNames',var);
fn0=fullfile(path_save,sprintf('%s-历史仓位.csv',key_str));
writetable(X,fn0)

for i = 1:length(re2)
    re2{i} = 1+re2{i};
end
sta_re = curve_static_batch(re2,index_info); %bog方法参数统计
sta_re2 = curve_static_batch(re6,cellfun(@(x) [x,'-LS'],index_info,'UniformOutput',false)); % bog方法对冲指数
re4=[re4{:}]'; %bog方法每年收益统计
re5=[re5{:}]'; %bog方法对冲指数每年收益统计

re7=[re7{:}]';
re7 = cell2table(re7,'VariableNames',{'tradeDate','Index','com_ticker','buyPrice','ma','lo'});
fn1 = fullfile(path_save,sprintf('%s-信号计算数值.csv',key_str));
writetable(re7,fn1);

