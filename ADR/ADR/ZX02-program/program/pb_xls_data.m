function sub_data=pb_xls_data(sub_fn)
    x=readtable(sub_fn);    
    x([1:2,4,6],:) = [];
    x=table2cell(x);
    x(1,:) = fillmissing(x(1,:),'previous');


    tref = x(3:end,1);
    tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));
    x=x(:,2:end);
    var_info = unique(x(2,:),'stable');

    stocks = x(1,:);
    info =x(2,:);
    x(1:2,:) = [];
    x = cellfun(@str2double,x);

    sub_data=[];
    s_re= cell(size(var_info));
    for i = 1:length(var_info)
        sub_ind = strcmp(info,var_info(i));
        sub_x = x(:,sub_ind);
        sub_data.(var_info{i}) = sub_x;
        s_re{i} = stocks(sub_ind);
    end
    sub_data.stocks=s_re{1};
    sub_data.tref = tref;
end