function conna = mysql_conn()
fid = fopen('para.json');
if fid>0
    p = fgetl(fid);
    fclose(fid);
    p=jsondecode(p);
    user_p = p.mysql_para.user_name;
    pass_wd = p.mysql_para.pass_wd;
    port = p.mysql_para.port;
else
    user_p = 'root';
    pass_wd = 'liudehua';
    port = 3306;
end
mysql_link = sprintf('jdbc:mysql://localhost:%d/futuredata?useSSL=false&',port);
conna = database('yuqerdata',user_p,pass_wd,'com.mysql.jdbc.Driver',mysql_link);