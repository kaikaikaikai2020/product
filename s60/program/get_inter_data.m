%不并行，直接使用in的方法
function x = get_inter_data(stocks,t0,tt)
    sql_str = ['select ticker,tradeDate,openPrice as op,highestPrice as hi,lowestPrice as lo,',...
    'closePrice as cl from yuqerdata.yq_MktEqudAdjAfGet where tradeDate>="%s"  and tradeDate<="%s" order by tradeDate'];    
    X = fetchmysql(sprintf(sql_str,t0,tt),2);
    ia = cellfun(@(x) any(ismember(stocks,x)),X(:,1));
    X = X(ia,:);

    tref = unique(X(:,2));
    tday = cellfun(@(x) str2double(strjoin(strsplit(x,'-'),'')),tref);
    symbol = unique(X(:,1));
    T = length(symbol);
    x=cell(T,1);
    
    for i = 1:T
        temp = X(strcmp(X(:,1),symbol(i)),2:end);
        [~,ia,ib] = intersect(tref,temp(:,1));
        sub_x = nan(length(tref),size(temp,2)-1);
        sub_x(ia,:) = cell2mat(temp(ib,2:end));
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
    x=[];
    x.op = op;
    x.hi=hi;
    x.lo=lo;
    x.cl=cl;
    x.stocks = symbol;
    x.tday = tday;
    x.tref = tref;
end