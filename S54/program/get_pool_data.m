function x = get_pool_data(stocks,t0,pool)
    if nargin <3
        pool = 'us';
    end
    if nargin <2
        t0 = '2006-05-11';
    end
    key_str = 'load index com data';
    if strcmp(pool,'us')
        sql_str = ['select tradeDate,openPrice as op,highPrice as hi,lowPrice as lo,',...
        'closePrice as cl from polygon.usastock_day where ticker="%s"and tradeDate>="%s" order by tradeDate'];
    elseif strcmp(pool,'US')
        %MktUsequdGetS54
        sql_str = ['select tradeDate,openPrice as op,highestPrice as hi,lowestPrice as lo,',...
        'closePrice as cl from yuqerdata.MktUsequdGetS54 where ticker="%s"and tradeDate>="%s" order by tradeDate'];
    elseif strcmpi(pool,'japan')
        sql_str = ['select tradeDate,openPrice as op,highestPrice as hi,lowestPrice as lo,',...
        'closePrice as cl from quandl.daytick_japan where ticker="%s"and tradeDate>="%s" order by tradeDate'];
    elseif strcmpi(pool,'HK')
        sql_str = ['select tradeDate,openPrice as op,highestPrice as hi,lowestPrice as lo,',...
        'closePrice as cl from yuqerdata.MktHKEqudGetS54 where ticker="%s"and tradeDate>="%s" order by tradeDate'];
    else
        sql_str = ['select tradeDate,openPrice as op,highestPrice as hi,lowestPrice as lo,',...
        'closePrice as cl from yuqerdata.yq_dayprice where symbol="%s"and tradeDate>="%s" order by tradeDate'];
    end
    T = length(stocks);
    X = cell(T,1);
    ind=true(size(stocks));
    parfor i = 1:length(stocks)
        sub_x1 = fetchmysql(sprintf(sql_str,stocks{i},t0),2);
        if isempty(sub_x1)
            ind(i) = false;
        end
        X{i} = sub_x1;
        if eq(i/100,floor(i/100))
            sprintf('%s loading data %d-%d',key_str,i,T)
        end
    end
    stocks = stocks(ind);
    X(cellfun(@isempty,X)) = [];
    T = length(X);
    tref = cellfun(@(x) x(:,1)',X,'UniformOutput',false);
    tref = unique([tref{:}]);
    tday = cellfun(@(x) str2double(strjoin(strsplit(x,'-'),'')),tref);
    x=cell(T,1);
    for i = 1:T
        temp = X{i};
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
    x=[];
    x.op = op;
    x.hi=hi;
    x.lo=lo;
    x.cl=cl;
    x.stocks = stocks;
    x.tday = tday;
    x.tref = tref;
end