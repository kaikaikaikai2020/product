%国内成分计算
clear;
%难点 S&P 500成分股代号获取
addpath(genpath(fullfile(pwd,'jplv7')))
topN=10; % Max number of positions
entryZscore=1;
lookback=20; % for MA

key_str = 'A_41';
index_pool = {'000300','000905','000001'};
index_info = {'I000300','I000905','I000001'};
T = length(index_pool);
re = cell(T,1);
tref0 = [];
for index_id = 1:3
    index_sel = index_pool{index_id};
    stocks = yq_methods.get_index_pool(index_sel,'2005-01-01');

    temp = get_pool_data(stocks,'2005-01-01','csi');
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

    retO2C=(cl-op)./op;


    pnl=smartsum(positionTable.*(retO2C), 2);
    ret=pnl/topN; 
    ret(isnan(ret))=0;
    
    
    sub_pool = cell(size(positionTable(:,1)));
    for i = 1:length(sub_pool)
        temp = positionTable(i,:);
        sub_sym = stocks(temp>0);
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
        title(sprintf('A-41-%s',index_sel))
    else
        figure;
        plot(cumprod(1+ret)-1); % Cumulative compounded return
    end
    
    pos0=[tref,sub_pool];
    pos0=cell2table(pos0,'VariableNames',{'date',index_info{index_id}});
    tref0=unique([tref0;table2cell(pos0(:,1))]);
    re{index_id}=pos0;
end


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

writetable(X,'A_41_bog_csi.csv')