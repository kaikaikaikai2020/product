%不并行，直接使用in的方法
function x = get_pool_data_update2(stocks,t0,pool)
    if nargin <3
        pool = 'us';
    end
    if nargin <2
        t0 = '2006-05-11';
    end
    key_str = 'load index com data';
    if strcmpi(pool,'us')
        sql_str = ['select tradeDate,openPrice as op,highPrice as hi,lowPrice as lo,',...
        'closePrice as cl from polygon.usastock_day where ticker in(%s) and tradeDate>="%s" order by tradeDate'];
    else
        sql_str = ['select tradeDate,openPrice as op,highestPrice as hi,lowestPrice as lo,',...
        'closePrice as cl from yuqerdata.yq_dayprice where symbol in (%s) and tradeDate>="%s" order by tradeDate'];
    end
    symbols = sprintf('"%s"',strjoin(stocks,'","'));
    
    
    X = fetchmysql(sprintf(sql_str,symbols),2);
    
    T = length(stocks);
    X = cell(T,1);
    parfor i = 1:length(stocks)
        X{i} = fetchmysql(sprintf(sql_str,stocks{i},t0),2);
        sprintf('%s loading data %d-%d',key_str,i,T)
    end
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