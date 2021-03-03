%clear
function  updateA38data()
    sprintf('begin update bloomberg data succeed')
    fn='bloomberg updates.xlsx';
    x=readtable(fn);
    x=x([3,7:end],:);
    ind = table2cell(x(2,:));
    ind = cellfun(@isnumeric,ind);
    x = x(:,~ind);
    data = [];
    data=formatDataS54(x(:,4:14),data);
    data=formatDataS54_f2(x(:,15:24),data);
    data=formatDataS54(x(:,25:31),data);
    data=formatDataS54(x(:,1:3),data);
    save('bloombergdata.mat','data')
    sprintf('%s',strjoin(fieldnames(data),','))
    sprintf('update bloomberg data succeed')
end

