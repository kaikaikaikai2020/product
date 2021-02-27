%{
支线任务2数据对接程序，请务必保证
数据的格式不要发生变化!
%}

function data = updatedataZX02()
    [~,~,X] = xlsread('Non Fungible ADR strategies.xlsx');
    c_p = X(4,:);
    c_p = c_p(cellfun(@isnumeric,c_p));
    c_p = cell2mat(c_p);
    c_p = c_p(~isnan(c_p));
    X = X([8,10,12:end],:);
    ind =cellfun(@isnumeric, X(2,:));
    X(:,ind) = [];

    ind = strcmpi(X(2,:),{'Premium (buy US close and sell Asia Open)'});
    X(:,ind) = [];
    ind = strcmpi(X(2,:),{'Premium (sell Asia close and buy US open)'});
    X(:,ind) = [];

    [~,n] = size(X);
    T = n/16;

    data = cell(T,1);

    for i = 1:T
        n1 = 16*(i-1)+1;
        n2 = 16*i;

        x = X(:,n1:n2-1);

        var_info = x(2,2:6);
        index = x(1,2:end);
        del_ind = cellfun(@isnumeric,index);
        index(del_ind) = [];
        index = cellfun(@(x) strsplit(x,' '),index,'UniformOutput',false);
        p_v = cellfun(@(x) x{2},index,'UniformOutput',false);
        index = cellfun(@(x) x{1},index,'UniformOutput',false);
        index = cellfun(@(x) replace(x,'+','_'),index,'UniformOutput',false);
        tref = x(3:end,1);

        %数字日期转换/以防万一
        ind = cellfun(@isnumeric,tref);
        if any(ind)
            tref(ind) = cellstr(datestr(cell2mat(tref(ind))+693960,'yyyy-mm-dd'));
        end


        tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));
        x= cell2mat(x(3:end,2:end));
        tempdata = [];
        for j = 1:3
            sub_ind1 = (j-1)*5+1;
            sub_ind2 = j*5;
            if sub_ind2 > size(x,2)
                sub_ind2 = size(x,2);
            end

            sub_x= x(:,sub_ind1:sub_ind2);
            temp = [];
            temp.data = sub_x;
            temp.tref = tref;
            tempdata.(sprintf('x%d',j)) = temp;

        end
        tempdata.index = index;
        tempdata.var_info = var_info;
        tempdata.p_v=p_v;
        tempdata.c_p = c_p(i);
        data{i} = tempdata;
    end
end