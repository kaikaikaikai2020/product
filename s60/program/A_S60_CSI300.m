clear
%% 生成股票对
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
%% 计算曲线
dos('python M_S60_indexTest.py')

%% 合成曲线
pn='S60P3_para_CSI300';
fns = dir(fullfile(pn,'*.csv'));
fns = {fns.name}';
T = length(fns);
tref0= yq_methods.get_tradingdate();
tref0_num = datenum(tref0);

t_all = cell(T,1);
for i = 1:T
    sub_title =split(fns{i},'.');
    sub_title = sub_title{1};
    sub_num = strsplit(sub_title,'-');
    t_all{i} = sub_title(end-9:end);
end

t_all = sort(unique(t_all));
T1 = length(t_all);
re = cell(T1,1);
for i = 1:T1
    sub_t1 = t_all{i};
    if i <T1
        sub_tt = t_all{i+1};
    else
        sub_tt = datestr(now,'yyyy-mm-dd');
    end
    ind1 = strcmp(t_all,sub_t1);
    %文件名
    sub_fns = fns(ind1);
    sub_tref0 = tref0(tref0_num>datenum(sub_t1) & tref0_num<=datenum(sub_tt));
    sub_re = zeros(length(sub_tref0),length(sub_fns));
    for j = 1:length(sub_fns)
        [~,~,x]= xlsread(fullfile(pn,sub_fns{j}));
        sub_tref = x(2:end,2);
        sub_tref = cellstr(datestr(datenum(sub_tref),'yyyy-mm-dd'));
        yc = cell2mat(x(2:end,3:4));
        yc = yc(:,1)-yc(:,2);
        
        [~,ia,ib] = intersect(sub_tref0,sub_tref);
        sub_re(ia,j) = yc(ib);
    end
    sub_re = [sub_tref0,num2cell(mean(sub_re,2))];
    re{i} = sub_re'; 
    sprintf('loading results %d-%d',i,T1)
end
Y = [re{:}]';
tref = Y(:,1);
yc = cumsum(cell2mat(Y(:,2)));
tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));
h = figure_S53(yc,tref,'000300组合');
