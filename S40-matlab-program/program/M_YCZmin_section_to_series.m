%{
��Ԥ���߽�������ת��Ϊʱ����������
���ڹ�������ݵ����ݿ⣬����Խ��Խ�󣬱������Խ��Խ��
��Ҫ���������ܹ��ж�ʱ���
%}

clear
keystr = 'Ԥ���� ����תʱ������';
dN= 'ycz_min_series';

dN_source = 'ycz_min_history';
symbols = yq_methods.get_symbol_A();
ycz_symbols = symbols;
id1 = cellfun(@(x) strcmp(x(1),'6'),symbols);
ycz_symbols(id1) = cellfun(@(x) ['sh',x],ycz_symbols(id1),'UniformOutput',false);
ycz_symbols(~id1) = cellfun(@(x) ['sz',x],ycz_symbols(~id1),'UniformOutput',false);

%��ȡԤ����section����
tns_sources = fetchmysql(sprintf('show tables from %s',dN_source),2);
tns_sources = sort(tns_sources);

T_tns_source = length(tns_sources);
sql_str = 'select distinct(symbol) from %s.`%s` order by symbol';
sql_str_f1 = 'insert into %s.%s select * from %s.`%s`  where symbol = ''%s''';

if exist('S40_Astock_min_date.mat','file')
    load('S40_Astock_min_date');
    num0 = find(strcmp(tns_sources,sub_tns_sources))+1;
else
    sub_tns_sources = [];
    num0 = 1;
end
for i = num0:T_tns_source
    sub_tns_sources = tns_sources{i};
    sub_symbols = fetchmysql(sprintf(sql_str,dN_source,sub_tns_sources),2);
    T_symbols = length(sub_symbols);
    parfor j = 1:T_symbols
        %����mysql���
        sub_sql_str = sprintf(sql_str_f1,dN,sub_symbols{j},dN_source,sub_tns_sources,sub_symbols{j});
        try
            exemysql(sub_sql_str);
        catch e_info
            sprintf(e_info.message)
        end
        sprintf('%s : %d-%d  %d -%d',keystr,j,T_symbols,i,T_tns_source)
    end     
end
if ~isempty(sub_tns_sources)
    save('S40_Astock_min_date.mat','sub_tns_sources');
end
