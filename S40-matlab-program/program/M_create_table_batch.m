%{
Ҫ��
1�������ǵ����ݿ��csi300������������
2����ͬ�����߼�����һ��ȫa�ɵ�1min�ı��֣������ǵ����ݣ�
3����ͬ�����߼�����һ��������ݣ���������ֵ��鷳����ÿ���1����Ӧ����1440�����ݣ�
                          �����õ�������ÿ�����ݵĸ�����̫һ��������Ҫ����һ��
4����ͬ�����߼�����һ��������Ʒ�ڻ�������
5��������иĽ��ĵط�����ӭ���Ľ�

Ϊ�˱��ڿ���ʱ���������ݣ������к�����A������ת��ʱ����������
%}
clear

%��ȡA������
%�������

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