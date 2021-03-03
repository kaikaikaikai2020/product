function [x,var_info]=trans_S54_data(x)
    info = containers.Map({'日期','收盘','开盘','高','低','交易量','涨跌幅'},...
        {'tradeDate','closePrice','openPrice','highestPrice','lowestPrice','volume','CHG'});

    var_info = cell(size(x(1,:)));
    for i = 1:size(x,2)
        var_info{i} = info(x{1,i});
    end

    x= x(2:end,:);
    x(:,1) = cellfun(@(x) strrep(x,'年','/'),x(:,1),'UniformOutput',false);
    x(:,1) = cellfun(@(x) strrep(x,'月','/'),x(:,1),'UniformOutput',false);
    x(:,1) = cellfun(@(x) strrep(x,'日',''),x(:,1),'UniformOutput',false);

    rm_s = {',','%'};
    for i = 2:size(x,2)
        for j = 1:length(rm_s)
            x(:,i) = cellfun(@(x) strrep(x,rm_s{j},''),x(:,i),'UniformOutput',false);
        end
    end

    ind1 = find(strcmp(var_info,'tradeDate'));
    ind2 = find(strcmp(var_info,'volume'));

    ind3 = 1:length(var_info);
    ind3([ind1,ind2]) = [];
    x(:,ind3) = cellfun(@str2double,x(:,ind3),'UniformOutput',false);
    % for i = 1:size(x,1)
    %     sub_str = x{i,ind2};
    %     if  contains(sub_str,'K') ||  contains(sub_str,'k')
    %         sub_str = strrep(sub_str,'K','');
    %         sub_str = strrep(sub_str,'k','');
    %         sub_str = str2double(sub_str)*1000;
    %     elseif contains(sub_str,'M') ||  contains(sub_str,'m')
    %         sub_str = strrep(sub_str,'K','');
    %         sub_str = strrep(sub_str,'k','');
    %         sub_str = str2double(sub_str)*1000;
    %     else
    %         sub_str = str2double(sub_str);
    %     end
    %     x{i,ind2} = sub_str;        
    % end
    tref_num = datenum(x(:,1));
    x(:,1) = cellstr(datestr(tref_num,'yyyy-mm-dd'));
end