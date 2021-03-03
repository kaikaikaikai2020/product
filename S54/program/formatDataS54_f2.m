function data = formatDataS54_f2(x1,data)
    if nargin <2
        data = [];
    end
    if istable(x1)
        x1 = table2cell(x1);
    end
    
    %x1 = x1([1,5:end],:);
    var1 = x1(1,2:2:end);
    var1 = cellfun(@(x) strrep(x,' Curncy',''),var1,'UniformOutput',false);
    % Index
    var1 = cellfun(@(x) strrep(x,' Index',''),var1,'UniformOutput',false);
    var1 = cellfun(@(x) strrep(x,'+',''),var1,'UniformOutput',false);
    x1 = x1(2:end,:);
    
    for i = 1:length(var1)
        temp = x1(:,i*2-1:i*2);
        del_ind = cellfun(@isempty,temp(:,1));
        temp = temp(~del_ind,:);
        temp(:,1) = cellstr(datestr(datenum(temp(:,1)),'yyyy-mm-dd'));
        temp(:,2) = cellfun(@str2double, temp(:,2:end),'UniformOutput',false);
        del_ind = cellfun(@iscell,temp(:,2));
        temp(del_ind,:) = [];
        data.(var1{i}) =temp;
    end
end