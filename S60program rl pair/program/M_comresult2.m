%��������
%ǰһ�����ݼ����������һ�����ݻز�
clear
%pn = 'HA_result';
%pn='AIndex_result';
%pn = 'AD_result';
pn='S60P3_para';
fns = dir(fullfile(pn,'*.csv'));
fns = {fns.name}';
T = length(fns);

obj = wordcom(fullfile(pwd,sprintf('S60%s.doc',pn)));

for i = 1:T
    sub_title =split(fns{i},'.');
    sub_title = sub_title{1};
    sub_fn = fullfile(pn,fns{i});
    
    [~,~,x]= xlsread(sub_fn);
    tref = x(2:end,2);
    yc = cell2mat(x(2:end,3:4));
    yc = yc(:,1)-yc(:,2);
    tmp = ceil(length(tref)*0.5);
    tref = tref(tmp:end);
    yc = yc(tmp:end);
    
    yc = cumsum(yc);
    tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));
    
    h = figure_S53(yc,tref,sub_title);
    pasteFigure(obj,h,sub_title);
    pause(1)
end

 CloseWord(obj)
