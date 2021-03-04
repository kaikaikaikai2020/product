%{
要求
1）用我们的数据库的csi300的数据来测试
2）用同样的逻辑测试一下全a股的1min的表现（用我们的数据）
3）用同样的逻辑测试一下外汇数据（外汇数据又点麻烦，它每天的1分钟应该是1440个数据，
                          但我拿到的数据每天数据的个数不太一样，所以要处理一下
4）用同样的逻辑测试一下其他商品期货的数据
5）如果您有改进的地方，欢迎您改进

为了便于控制时间序列数据，将现有横截面的A股数据转换时间序列数据
%}
clear

%获取A股数据
%创建表格

symbols = yq_methods.get_symbol_A();
ycz_symbols = symbols;
id1 = cellfun(@(x) strcmp(x(1),'6'),symbols);
ycz_symbols(id1) = cellfun(@(x) ['sh',x],ycz_symbols(id1),'UniformOutput',false);
ycz_symbols(~id1) = cellfun(@(x) ['sz',x],ycz_symbols(~id1),'UniformOutput',false);
dN= 'ycz_min_series';

for i = 1:length(ycz_symbols)
    tn = ycz_symbols{i};
    var_info={'symbol', 'tradingdate', 'open', 'high', 'low', 'close', 'volume', 'amount','volume2'};
    %var_info=['varchar(8)','datetime','float','float','float','float','float','float','float']
    %var_info = {'tradingdate', 'method_ID', 'index_code', 'more_r', 'less_r'};
    var_type = cell(size(var_info));
    var_type(:) = {'float'};
    var_type(1:2) = {'varchar(8)','datetime'};
    %key_var = {'symbol','tradingdate'};
    key_var = strjoin(var_info([1,2]),',');
    %key_var = var_info{1};
    create_table_adair(dN,tn,var_info,var_type,key_var)
end