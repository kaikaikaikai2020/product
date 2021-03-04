%方法修正
%前一半数据计算参数，后一半数据回测
clear
%pn = 'HA_result';
%pn='AIndex_result';
%pn = 'AD_result';
pn='S60P3_para';
fns = dir(fullfile(pn,'*.csv'));
fns = {fns.name}';
[~,ia] = intersect(fns,{'bac-ccfx_ic-ccfx_if.csv','bac-SPY-QQQ.csv','bac-ES_day-VX_day.csv',...
    'bac-ES_day-NQ_day.csv','bac-ES-NQ.csv',...
    'bac-ES-VX.csv',''});
fns = fns(ia);
T = length(fns);

obj = wordcom(fullfile(pwd,sprintf('S60%s_SPYQQQ.doc',pn)));

for i = 1:T
    sub_title =split(fns{i},'.');
    sub_title = sub_title{1};
    sub_num = strsplit(sub_title,'-');
    sub_num = str2double(sub_num{end});
    
    sub_fn = fullfile(pn,fns{i});
    
    [~,~,x]= xlsread(sub_fn);
    tref = x(2:end,2);
    yc = cell2mat(x(2:end,3:4));
    yc = yc(:,1)-yc(:,2);
    if isnan(sub_num)
        sub_num = floor(length(tref)/2);
    end
    tref = tref(sub_num:end);
    yc = yc(sub_num:end);
    
    yc = cumsum(yc);
    tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));
    
    h = figure_S53(yc,tref,sub_title);
    pasteFigure(obj,h,sub_title);
    pause(1)
end

 CloseWord(obj)
