function x = list_datelimit(t_str,window1)
if nargin < 2
    window1 = 180;
end

sql_str4 = ['select ticker,listDate from yuqerdata.equget where listStatusCd !=''UN''',...
                    'and listDate is not null'];
symbol_info = fetchmysql(sql_str4,2);
symbol_listdate = datenum(symbol_info(:,2));
ind = datenum(t_str)-symbol_listdate>window1;
x = symbol_info(ind,1);