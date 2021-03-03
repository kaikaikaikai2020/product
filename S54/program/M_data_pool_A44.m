%{
1）得到过去10天的开盘竞价移动平均
2)  得到开盘竞价的成交量（volume）
3）过滤出最大成交金额（开盘竞价移动平均*成交量）的N只股票
4）用这N只股票形成股票池，用Andrewlo part2 算法计算收益
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
