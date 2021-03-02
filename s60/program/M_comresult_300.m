%方法修正
%前一半数据计算参数，后一半数据回测
clear
%pn = 'HA_result';
%pn='AIndex_result';
%pn = 'AD_result';
pn='S60P3_300';
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
    
end
Y = [re{:}]';
tref = Y(:,1);
yc = cumsum(cell2mat(Y(:,2)));
tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));
h = figure_S53(yc,tref,'000300组合');
