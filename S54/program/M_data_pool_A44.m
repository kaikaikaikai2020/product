%{
1���õ���ȥ10��Ŀ��̾����ƶ�ƽ��
2)  �õ����̾��۵ĳɽ�����volume��
3�����˳����ɽ������̾����ƶ�ƽ��*�ɽ�������Nֻ��Ʊ
4������Nֻ��Ʊ�γɹ�Ʊ�أ���Andrewlo part2 �㷨��������
%}

fn='HSI open auction volume.xlsx';
[~,~,x] = xlsread(fn);
V = x(867:end,2:end);
tref = x(867:end,1);
symbol = x(4,2:end);
V(~cellfun(@isnumeric,V)) = {nan};
V=cell2mat(V);
tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));

symbol = cellfun(@(x) del_str(x),symbol,'UniformOutput',false);


data = [];
data.HSI.V = V;
data.HSI.tref = tref;
data.HSI.symbol = symbol;

fn='HSCEI open auction volume.xlsx';
[~,~,x] = xlsread(fn);
V = x(873:end,2:end);
tref = x(873:end,1);
symbol = x(4,2:end);
V(~cellfun(@isnumeric,V)) = {nan};
V=cell2mat(V);

tref = cellstr(datestr(datenum(tref),'yyyy-mm-dd'));
symbol = cellfun(@(x) del_str(x),symbol,'UniformOutput',false);

data.HSCEI.V = V;
data.HSCEI.tref = tref;
data.HSCEI.symbol = symbol;

save HK_vol data

function x = del_str(x)
x =strsplit(x,' ');
x = x{1};
x = sprintf('%0.5d',str2double(x));
end
