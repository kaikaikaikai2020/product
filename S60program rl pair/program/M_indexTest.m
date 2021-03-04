clear
%时间点选择
num_pair=5;
sql_str = ['select endDate from yuqerdata.yq_index_month ',...
    'where symbol = "000300" and endDate>="2015-01-01" order by endDate'];
index = '000300';
tref=fetchmysql(sql_str,2);
%在每个时间点执行筛选
T = length(tref);
para_pool = cell(T,1);
parfor i = 2:T
    sub_t0=tref{i-1};
    sub_t = tref{i};
    %com
    ticker = yq_methods.get_index_pool(index,sub_t);
    %time st
    ticker_st = yq_methods.get_stpt_symbol(sub_t);
    ticker = setdiff(ticker,ticker_st);
    %time
    ticker_dt = yq_methods.get_time_cut_symbol(datestr(datenum(sub_t)-10*365,'yyyy-mm-dd'));    
    ticker = intersect(ticker,ticker_dt);
    
    x = get_inter_data(ticker,sub_t0,sub_t);
    ticker = x.stocks;
    x = x.cl;
    ind = isnan(sum(x));
    ticker(ind) = [];
    x(:,ind) = [];
    r = corr(x);
    r = tril(r);
    r(r==1) = 0;
    rv=sort(r(:),'descend');
    sub_re = cell(num_pair,1);
    for j = 1:num_pair
        [ia,ib] = find(r==rv(j));
        temp = cell(size(ia));
        for k = 1:length(ia)
            temp{k} = [ia(k),ib(k)]';
        end
        sub_re{j} = [temp{:}];
    end
    sub_re = [sub_re{:}]';
    sub_re = sub_re(1:num_pair,:);
    
    sub_pool = ticker(sub_re(:,[1:end,end]));
    sub_pool(:,end) = {sub_t};
    para_pool{i} = sub_pool';
end

para_pool = [para_pool{:}]';
%save para_pool000300 para_pool
para_pool = cell2table(para_pool,'VariableNames',{'ticker1','ticker2','t'});
writetable(para_pool,'T000300.csv');