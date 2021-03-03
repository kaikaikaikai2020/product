function data = formatDataS54(x1,data)
    if nargin <2
        data = [];
    end
    if istable(x1)
        x1 = table2cell(x1);
    end
    %x1 = x1([1,5:end],:);
    var1 = x1(1,2:end);
    var1 = cellfun(@(x) strrep(x,' Curncy',''),var1,'UniformOutput',false);
    % Index
    var1 = cellfun(@(x) strrep(x,' Index',''),var1,'UniformOutput',false);
    var1 = cellfun(@(x) strrep(x,'+',''),var1,'UniformOutput',false);
    x1 = x1(2:end,:);
    del_ind = cellfun(@isempty,x1(:,1));
    x1 = x1(~del_ind,:);
    tref = cellstr(datestr(datenum(x1(:,1)),'yyyy-mm-dd'));
    x1 = cellfun(@str2double, x1(:,2:end));
    for i = 1:length(var1)
        temp = [tref,num2cell( x1(:,i))];
        del_ind = isnan(x1(:,i));
        temp(del_ind,:) = [];
        data.(var1{i}) =temp;
    end
end