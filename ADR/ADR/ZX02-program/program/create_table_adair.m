function [OK1,OK2,OK3] = create_table_adair(dN,tn,var_info,var_type,key_var)
    dN_all = fetchmysql('show databases;',2);
    if istable(dN_all)
        dN_all = table2cell(dN_all);
    end
    if ~any(strcmpi(dN_all,dN))
        exemysql(sprintf('create database %s',dN));
    end
    %check tables `
    tns_all = fetchmysql(sprintf('show tables from %s',dN),2);
    if istable(tns_all)
        tns_all = table2cell(tns_all);
    end

    if ~any(strcmpi(tns_all,tn))
        %create table   
        obj = mysqlTool();
        sqlquery1=obj.createTable(dN,tn,var_info,var_type);
        OK1 = exemysql(sqlquery1);
        OK2 = exemysql(sprintf('alter table %s.%s engine=MyISAM;',dN,tn));
        if ~isempty(key_var)
            OK3 = exemysql(sprintf('alter table %s.%s add primary key(%s);',dN,tn,key_var));
        else
            OK3 = true;
        end
    else
        OK1 = true;
        OK2 = OK1;
        OK3 = OK1;
    end
end