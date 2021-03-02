clear
%��������
%�����źŲ����� ����������ļ���S60P3_para��
%����ļ�����������
    %csv �ز�����
    %excel �����ź�
    %pkl �����м����
dos('python M_S60_ON.py')

%
pn='S60P3_para';
fns = dir(fullfile(pn,'*-800.csv'));
fns = {fns.name}';
T = length(fns);

obj = wordcom(fullfile(pwd,sprintf('S60%s-%s.doc',pn,datestr(now,'yyyy-mm-dd'))));

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
    tref = tref(sub_num:end);
    yc = yc(sub_num:end);
    
    yc = cumsum(yc);
    tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));
    
    h = figure_S53(yc,tref,sub_title);
    pasteFigure(obj,h,sub_title);
    pause(1)
end

 CloseWord(obj)