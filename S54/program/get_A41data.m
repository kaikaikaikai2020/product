function stks = get_A41data()
    sql_str = 'select * from S54.A41data order by tradeDate';
    X = fetchmysql(sql_str,2);
    %tref = cellfun(@(x) x(:,1)',X,'UniformOutput',false);
    tref = unique(X(:,2));
    tday = cellfun(@(x) str2double(strjoin(strsplit(x,'-'),'')),tref);
    stocks = unique(X(:,1));
    T = length(stocks);
    x=cell(T,1);
    parfor i = 1:T
        ind = strcmp(X(:,1),stocks(i));
        temp = X(ind,2:end);
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
    stks=[];
    stks.op=temp{1};
    stks.hi=temp{2};
    stks.lo=temp{3};
    stks.cl=temp{4};
    stks.stocks= stocks;
    stks.tday=tday;
    stks.tref = tref;
end