%{
创建表格
%}
clear

dN= 'S40';
tn = 'A_stock_signal';
var_info={'symbol', 'tradenum', 'f_val'};
var_type(1:3) = {'varchar(8)','int','int'};
key_var = strjoin(var_info([1,2]),',');
%key_var = var_info{1};
create_table_adair(dN,tn,var_info,var_type,key_var)