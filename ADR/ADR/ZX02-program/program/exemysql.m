function OK = exemysql(sqlquery,dN,usr,pass)
if nargin < 2
    dN = 'futuredata';
end
if nargin < 3
    usr = 'root';
end
if nargin < 4
    pass = 'liudehua';
end
conna = mysql_conn();
try
    if iscell(sqlquery)
        for i = 1:length(sqlquery)
            exec(conna,sqlquery{i});
        end
    else
        exec(conna,sqlquery)
    end
    OK = true;
catch
    close(conna);
    OK = false;
end
close(conna);
