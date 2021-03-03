%国内成分计算
clear;
%难点 S&P 500成分股代号获取
addpath(genpath(fullfile(pwd,'jplv7')))
topN=10; % Max number of positions
entryZscore=1;
lookback=20; % for MA

key_str = 'A_41_japan';

load('index_japan');
load('bloombergdata.mat');
fee1 = 1/10000;
fee2 = 1/10000;

index_pool = index_info;
%index_info = {'I000300','I000905','I000001'};
T = length(index_pool);
re = cell(T,1);
tref0 = [];
re2 =cell(T,1);
re3 = cell(size(re2));
re4 = cell(T,1);
re5 = cell(T,1);
re6 = cell(T,1);
for index_id = 1:T
    stocks = index_ticker{index_id};

    temp = get_pool_data(stocks,'2012-01-01','japan');
    stocks = temp.stocks;
    cl=temp.cl;
    hi=temp.hi;
    lo=temp.lo;
    op=temp.op;
    tday=temp.tday;
    tref = temp.tref;


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

    %retO2C=(cl-op)./op;
    retO2C=(cl*(1-fee1)-op*(1+fee2))./(op*(1+fee2));

    pnl=smartsum(positionTable.*(retO2C), 2);
    %fee2 = smartsum(positionTable.*(fee*ones(size(retO2C))),2);
    
    ret=pnl/topN; 
    ret(isnan(ret))=0;
    
    %fee2 = fee2/topN;    
    %ret = ret-fee2;
    
    sub_pool = cell(size(positionTable(:,1)));
    sub_pool_num = zeros(size(sub_pool));
    for i = 1:length(sub_pool)
        temp = positionTable(i,:);
        sub_sym = stocks(temp>0);
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
    % APR=8.7%, Sharpe=1.5

    %cumret=cumprod(1+ret)-1; % compounded ROE
    if size(tref,2)>1
        tref = tref';
    end
    %plot(cumret);
    if ~isempty(tref)
        y_re = cumprod(1+ret)-1;
        %setfigure
        h = figure_S53(y_re,tref,[]);
        title(sprintf('A-41-%s',index_info{index_id}))
    else
        y_re = cumprod(1+ret)-1;
        figure;
        plot(y_re); % Cumulative compounded return
    end
    
    pos0=[tref,sub_pool];
    pos0=cell2table(pos0,'VariableNames',{'date',index_info{index_id}});
    tref0=unique([tref0;table2cell(pos0(:,1))]);
    re{index_id}=pos0;
    re2{index_id} = y_re;
    %{
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
            temp(:,1) = index_info(index_id);
            info2{i} = temp';
        end
    end
    info2 = [info2{:}];
    re3{index_id} = info2;
    %}
    temp = get_year_return(tref,ret);
    temp(:,1) = index_pool(index_id);
    
    re4{index_id} = temp';
    
    index=data.(index_pool{index_id});
    [tref2,ia,ib] = intersect(tref,index(:,1));
    index = cell2mat(index(ib,2));
    index(2:end) = index(2:end)./index(1:end-1)-1;
    index(1) = 0;
    
    ret=ret(ia,:);
    ret(1) = 0;
    sub_pool_num = sub_pool_num(ia);
    ret_ls=(ret-index.*sub_pool_num)./(1+sub_pool_num);
    y_ls=cumprod(1+ret_ls);    
    %对冲结果
    
    temp = get_year_return(tref2,ret_ls);
    temp(:,1) = {sprintf('%sls',index_pool{index_id})};    
    re5{index_id} = temp';
    re6{index_id} = y_ls;
    
    h = figure_S53(y_ls-1,tref2,sprintf('A41bog-%s-LS',index_pool{index_id}));
    
    
end

%{
re3 = [re3{:}]';
re3 = cell2table(re3,'VariableNames',{'code','date','code2','buy','sel'});
writetable(re3,'A_41_bog_japan_detail.csv')
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

writetable(X,'A_41_bog_japan.csv')

for i = 1:length(re2)
    re2{i} = 1+re2{i};
end
sta_re = curve_static_batch(re2,index_info);
sta_re2 = curve_static_batch(re6,cellfun(@(x) [x,'-LS'],index_info,'UniformOutput',false));
re4=[re4{:}]';
re5=[re5{:}]';