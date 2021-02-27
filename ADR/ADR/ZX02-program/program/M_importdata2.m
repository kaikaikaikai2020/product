clear
[~,~,X] = xlsread('data2.xlsx');
[m,n] = size(X);
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
    data{i} = tempdata;
    
end
save data20201030 data
