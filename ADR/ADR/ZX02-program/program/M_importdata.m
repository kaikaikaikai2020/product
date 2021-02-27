clear
[~,~,x] = xlsread('data.xlsx');
data=[];
var_info = x(2,2:6);
index = unique(x(1,2:end),'stable');
index = cellfun(@(x) strsplit(x,' '),index,'UniformOutput',false);
index = cellfun(@(x) x{1},index,'UniformOutput',false);
index = cellfun(@(x) replace(x,'+','_'),index,'UniformOutput',false);

tref = x(3:end,1);
tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));
x= cell2mat(x(3:end,2:end));

for i = 1:3
    sub_ind1 = (i-1)*5+1;
    sub_ind2 = i*5;
    if sub_ind2 > size(x,2)
        sub_ind2 = size(x,2);
    end
    
    sub_x= x(:,sub_ind1:sub_ind2);
    temp = [];
    temp.data = sub_x;
    temp.tref = tref;
    data.(index{i}) = temp;
end


[~,~,x] = xlsread('data.xlsx','sheet2');
index = unique(x(1,2:end),'stable');
index = cellfun(@(x) strsplit(x,' '),index,'UniformOutput',false);
index = cellfun(@(x) x{1},index,'UniformOutput',false);
index = cellfun(@(x) replace(x,'+','_'),index,'UniformOutput',false);

tref = x(3:end,1);
tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));
x= cell2mat(x(3:end,2:end));

for i = 1:3
    sub_ind1 = (i-1)*5+1;
    sub_ind2 = i*5;
    if sub_ind2 > size(x,2)
        sub_ind2 = size(x,2);
    end
    
    sub_x= x(:,sub_ind1:sub_ind2);
    temp = [];
    temp.data = sub_x;
    temp.tref = tref;
    data.(index{i}) = temp;
end
data.var_info = var_info;

save data20201027 data
