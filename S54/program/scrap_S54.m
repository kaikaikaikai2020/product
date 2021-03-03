clear
%load testdata.mat
x = getTableFromWeb_mod_adair('https://cn.investing.com/indices/eu-stoxx50-historical-data', 2);
[x,var_info] = trans_S54_data(x);

tN = 'S54.A71data';
t0 = fetchmysql(sprintf('select max(tradeDate) from %s',tN),2);
tref_num = datenum(x(:,1));
ind = tref_num>=datenum(t0);
if any(ind)
    sub_x = x(ind,:);
    %删除最后一组数据
    exemysql(sprintf('delete from %s where tradeDate>="%s"',tN,t0{1}));
    datainsert_adair(tN,var_info,sub_x);
end

