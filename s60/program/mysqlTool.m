%sql server tool
%Methods
%1执行语句 function Message = dosql(sqlquery)
%2复制表结果但不复制表内容function sqlquery=copyTable(tabOld,tabNew,varnames)
%3新建表语句function sqlquery=createTable(dbname,tablename,varnames,Type)
%4删除表语句 function sqlquery=dropTable(tN)
%5新建数据库function sqlquery=createDatabase(dbName)
%6删除数据库 function sqlquery = deleteDatabase(dbName)
%7添加主键function sqlquery = addKey(dbName,tN,keyN,varnames)
%8删除主键function sqlquery = delKey(dbName,tN,keyN)
%9插入列function sqlquery = addColumn(dbName,tN,varnames,varType)
%10导入csv function sqlquery = importCSV(dbName,tN,fn,num)
%11function sqlquery = removeDuplicate(dN,tN)
%12获取字段function sqlquery = getColumnName(dN,tN)
%13选择一张表的数据插入另一张表function sqlquery = insertContent(dN1,tN1,dN2,tN2)
%14清空表内容语句function sqlquery = truncateTable(dN,tN)
%15设置表格式cellarray or numeric sqlquery = setDataFormat(sel)
%16判断表格sqlquery = getTable(),1有，0无
%17将cell字段转换为字符，便于后续选择 sqlquery = transCellStr(varName)
%18 获取数据口下所有表 sqlquery = getTableName(dN)
classdef mysqlTool < handle
    methods(Static)
        %1执行语句
        function Message = dosql(sqlquery)
            conna = database('sql2008','','');
            try
                info = exec(conna,sqlquery);
                Message = info.Message;
                close(conna);
            catch
                Message = '程序出错';
                close(conna);
            end
        end
        %2复制表结果但不复制表内容
        function sqlquery=copyTable(tabOld,tabNew,varnames)
            if eq(nargin,3)
                str1 = [];
                for i = 1:length(varnames)
                    str1 = [str1,varnames{i},', '];
                end
                str1(end-1:end) = [];
            else
                str1 = '*';
            end
            sqlquery = ['select ',str1,' into ',tabNew,'  from ',tabOld,...
                ' where 1=2'];
        end
        %3新建表语句
        function sqlquery=createTable(dbname,tablename,varnames,Type)
            if eq(nargin,3)            	
                Type = cell(size(varnames));
                Type(1:end) = {'float'};
            end
            sqlquery = ['create table ',dbname,'.',tablename,'('];
            for i = 1:length(varnames)
                sqlquery = [sqlquery,' ',varnames{i},' ',Type{i},',',10];
            end
            sqlquery = [sqlquery(1:end-2),')'];
        end
        %4删除表语句
        function sqlquery=dropTable(tN)
            sqlquery = ['if exists(select * from sysobjects where name=''',tN,''')',10,...
            'drop table ',tN];
        end
        %5新建数据库
        function sqlquery=createDatabase(dbName)
            sqlquery = ['create database ',dbName,10,... 
            ' on  primary  -- 默认就属于primary文件组,可省略 ',10,...
            '(',10,...
            '/*--数据文件的具体描述--*/',10,...
                'name=''' dbName,'_data''',',  -- 主数据文件的逻辑名称',10,...
                'filename=','''D:\sqlData\', dbName,'_data.mdf''',10,...
            ') ',10,...
            'log on',10,...
            '(',10,...
            '/*--日志文件的具体描述,各参数含义同上--*/',10,...
                'name=''',dbName,'_log''',',',10,...
                'filename=''','D:\sqlData\',dbName,'_log.ldf''',10,...
            ')'];
        end
        %6删除数据库
        function sqlquery = deleteDatabase(dbName)
            sqlquery=['use master -- 设置当前数据库为master,以便访问sysdatabases表',10,...
            'go',10,...
            'if exists(select * from sysdatabases where name=''',dbName,''')',10,...
            'drop database ',dbName,10,...
            'go'];
        end
        %7添加主键
        function sqlquery = addKey(dbName,tN,keyN,varnames)
            sqlquery=['alter table ',dbName,'.',tN,' add constraint ',keyN,' primary key('];
            for i = 1:length(varnames)
                sqlquery = [sqlquery,varnames{i},', '];
            end
            sqlquery(end-1:end) = ' )'; 
            %alter table dbo.GData1 add constraint pk_s primary key (stkn,stktt)
        end
        %8删除主键
        function sqlquery = delKey(dbName,tN,keyN)
            %alter table 表名 drop constraint 主键名 sp_help 表
            sqlquery=['use ',dbName,10,...
                ' alter table ',tN, ' drop constraint ',keyN];
        end
        %9插入列
        function sqlquery = addColumn(dbName,tN,varnames,varType)
            sqlquery = ['use ',dbName,10,...
                ' alter table ',tN,' add '];
            if iscell(varnames)
                for i = 1:lenth(varnames)
                    sqlquery = [sqlquery, varnames{i} ,' ',varType{i},', ']; 
                end
                sqlquery(end-1:end) = [];
            else
                sqlquery = [sqlquery,varnames,' ',varType];
            end
            %Alter Table dbo.GData1 ADD var1 int null
        end
        %10导入csv
        function sqlquery = importCSV(dbName,tN,fn,num)
            if eq(nargin,3)
                num=1;
            end
            sqlquery = ['bulk insert ',dbName,'.dbo.',tN,' ',10,...
                'from ''',fn,'''',10,...
                'with(',10,...
                'fieldterminator = '','' ,',10,...
                'firstrow = ',num2str(num),',',10,...                
                'rowterminator=''0x0a'')'];
%              'ERRORFILE =''D:\sqlData\error.txt'',',10,...            
%             BULK INSERT CSVTable
%             FROM 'D:\csv.txt'
%             WITH(
%              FIELDTERMINATOR = ',',
%              ROWTERMINATOR = '\n'
%             )
%             SELECT * FROM CSVTable
        end
        %11表格去重
        function sqlquery = removeDuplicate(dN,tN)
            sqlquery = ['use ',dN,10,...
                ' select distinct * into #Tmp1 from ',tN,10,...
                ' drop table ',tN,10,...
                ' select * into ',tN,' from #Tmp1 ',10,...
                ' drop table #Tmp1 '];
            % select distinct * into #Tmp from tableName
            % drop table tableName
            % select * into tableName from #Tmp
            % drop table #Tmp
        end
        %12获取表字段
        function sqlquery = getColumnName(dN,tN)
            sqlquery = ['Select Name FROM SysColumns Where id=Object_Id(''',...
                dN,'.dbo.',tN,''')'];
        end
        %Select Name FROM SysColumns Where id=Object_Id('Tmp.dbo.stockData')
        %13选择一张表的数据插入另一张表
        function sqlquery = insertContent(dN1,tN1,dN2,tN2)
            sqlquery = ['insert into ',dN1,'.dbo.',tN1,10,...
                'select * from ',dN2,'.dbo.',tN2];
        end
        %insert into research.dbo.sec_priceV2
        %select * from Tmp.dbo.stockData
        %14清空表内容语句
        function sqlquery = truncateTable(dN,tN)
            sqlquery = ['truncate table ',dN,'.dbo.',tN];
        end
        %15设置表格式cellarray or numeric sqlquery = setDataFormat(sel)
        function sqlquery = setDataFormat(sel)
            if strcmp(sel,'cellarray')
                sqlquery='setdbprefs(''DataReturnFormat'',''cellarray'')';
                eval(sqlquery);
            elseif strcmp(sel,'numeric')
                sqlquery='setdbprefs(''DataReturnFormat'',''numeric'')';
                eval(sqlquery);
            else
                sqlquery = [];
                sprintf('Someting is wrong!')
                keyboard;
            end
        end
        %16判断表格sqlquery = getTable(),1有，0无
        %注意Tmp.dbo.infoTable必须已经存在，否则手动创建
        function sqlquery = getTable(dN,tN)
            sqlquery=['use ',dN,10,...
            'IF EXISTS  (SELECT  * FROM dbo.SysObjects WHERE ID = object_id(N''',...
            tN,''') AND OBJECTPROPERTY(ID, ''IsTable'') = 1)',10,...
            'select Y from Tmp.dbo.infoTable',10,...
            'ELSE',10,...
            'select N from Tmp.dbo.infoTable'];
        end
        %17将cell字段转换为字符，便于后续选择 sqlquery = transCellStr(varName)
        function sqlquery = transCellStr(varName)
            str1 = cellfun(@(x) [x,','],varName,'UniformOutput',false);
            str1 = [str1{:}];
            sqlquery = str1(1:end-1);
        end
        %18 获取数据口下所有表 sqlquery = getTableName(dN)
        function sqlquery = getTableName(dN)
             sqlquery = ['SELECT Name FROM ',dN,'..SysObjects Where XType=''U'' ORDER BY Name'];
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%获取主键名称
% SELECT A.NAME AS 表名,B.NAME AS 主键名 
% FROM  SYSOBJECTS A 
%     JOIN SYSOBJECTS B 
%         ON A.ID=B.PARENT_OBJ 
%         AND A.XTYPE='U' AND B.XTYPE='PK'
%剔除字段中回车符

%找重复
% use threeReportData
% update revenueV52
% SET stockID = REPLACE(stockID, CHAR(13), '')
%use Tmp
%select ticker from stockData
%where ticker in (
%select ticker from stockData
%group by ticker
%having(count(*))>1
%)
