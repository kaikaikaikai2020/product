%��yuqer��ȡ���ݷ����ϼ�
classdef yq_methods< handle
    methods(Static)
        function tref_m = get_month_data()
            sql_str = ['select endDate from yuqerdata.yq_index_month where ',...
                'symbol = ''000001'' order by endDate'];
            tref_m = fetchmysql(sql_str,2);
            
        end
        function symbol = get_symbol_A()
            sql_str = ['select distinct(ticker)  from yuqerdata.equget ',...
                'where equTypeCD = ''A'' and listStatusCD !=''UN'' and ',...
                'ListSectorCD<=3 and ticker not like ''D%'' order by ticker;'];
            symbol = fetchmysql(sql_str,2);
            
            
        end
        %����t1 ���  t2 �յ�
        function tref = get_tradingdate_future(t1,t2)
            if nargin < 1
                t1 = [];
            end
            if nargin < 2
                t2=[];
            end
            
            if isempty(t1)
                t1 = '1980-01-01';
            end
            if isempty(t2)
                t2 = datestr(now+200,'yyyy-mm-dd');
            end
            
            sql_str = ['select tradingdate from yuqerdata.yq_tradingdate_future ',...
                'where tradingdate >= ''%s'' and tradingdate<=''%s'' order by tradingdate ;'];
            sql_str = sprintf(sql_str,t1,t2);
            tref = fetchmysql(sql_str,2);
        end
        %
        function tref = get_tradingdate(t1,t2)
            if nargin < 1
                t1 = [];
            end
            if nargin < 2
                t2=[];
            end
            
            if isempty(t1)
                t1 = '1980-01-01';
            end
            if isempty(t2)
                t2 = datestr(now,'yyyy-mm-dd');
            end
            
            sql_str = ['select tradeDate from yuqerdata.yq_index ',...
                'where tradeDate >= ''%s'' and tradeDate<=''%s'' and ',...
                'symbol = ''000001'' order by tradeDate ;'];
            sql_str = sprintf(sql_str,t1,t2);
            tref = fetchmysql(sql_str,2);
        end
        %��ȡ�µ����һ��
        function [month_cut_date1,month_cut_date2] = get_month_day(tref)
            tref_num = datenum(tref);
            month_index = month(tref_num);
            month_cut = [0;find(diff(month_index))];
            month_cut = [month_cut(1:end-1)+1,month_cut(2:end)];
            month_cut_date1 = tref(month_cut(:,1));
            month_cut_date2 = tref(month_cut(:,2));
        end
        %��ȡҵ���챨����
        %symbol���������ڣ���ƽ������ڣ�ʱ�䳤�ȣ����ͣ���ֵ
        function x = get_YeJiKuaiBao(var_name)
            sql_str1 = ['select ticker,publishdate,enddate,fiscalPeriod,',...
                'reportType,%s from yuqerdata.yq_FdmtEeGet '];
            sql_str2 = [' where enddate>''2005-01-01'' ',...
                'and mergedFlag = 1 and (secID like ''0%'' or secID like ''60%'' or secID like ''30%'')  ',... 
                'and reporttype in (''Q1'',''S1'',''CQ3'',''A'') order by secID,enddate desc'];
            sql_str1 = sprintf(sql_str1,var_name);
            sql_str = [sql_str1,sql_str2];
            x = fetchmysql(sql_str,2);
        end
        %��ȡ�ϲ����������
        %symbol���������ڣ���ƽ������ڣ�ʱ�䳤�ȣ����ͣ���ֵ
        function x = get_HeBingLiRun(var_name)
            sql_str1 = ['select ticker,publishdate,enddate,fiscalPeriod,',...
                'reportType,%s from yuqerdata.nincome '];
            sql_str2 = [' where enddate>''2005-01-01'' and enddate=endDaterep ',...
                'and mergedFlag = 1 and (secID like ''0%'' or secID like ''60%'' or secID like ''30%'')  ',... 
                'and reporttype in (''Q1'',''S1'',''CQ3'',''A'') order by secID,enddate desc'];
            sql_str1 = sprintf(sql_str1,var_name);
            sql_str = [sql_str1,sql_str2];
            x = fetchmysql(sql_str,2);
        end
        %��ȡ�ϲ��ʲ���ծ��
        %symbol���������ڣ���ƽ������ڣ�ʱ�䳤�ȣ����ͣ���ֵ
        function x = get_HeBingZiChanFuZhai(var_name)
            sql_str1 = ['select ticker,publishdate,enddate,fiscalPeriod,',...
                'reportType,%s from yuqerdata.yq_FdmtBSGet '];
            sql_str2 = [' where enddate>''2005-01-01'' and enddate=endDaterep ',...
                'and mergedFlag = 1 and (secID like ''0%'' or secID like ''60%'' or secID like ''30%'')  ',... 
                'and reporttype in (''Q1'',''S1'',''Q3'',''A'') order by secID,enddate desc'];
            sql_str1 = sprintf(sql_str1,var_name);
            sql_str = [sql_str1,sql_str2];
            x = fetchmysql(sql_str,2);
        end
        %��ȡ�ϲ��ֽ�������
        %symbol���������ڣ���ƽ������ڣ�ʱ�䳤�ȣ����ͣ���ֵ
        function x = get_HeBingXianJinLiu(var_name)
            sql_str1 = ['select symbol,publishdate,enddate,fiscalPeriod,',...
                'reportType,%s from yuqerdata.yq_FdmtCFGetAll '];
            sql_str2 = [' where enddate>''2007-01-01'' and enddate=endDaterep ',...
                'and mergedFlag = 1 and (secID like ''0%'' or secID like ''60%'' or secID like ''30%'')  ',... 
                'and reporttype in (''Q1'',''S1'',''CQ3'',''A'') order by secID,enddate desc'];
            sql_str1 = sprintf(sql_str1,var_name);
            sql_str = [sql_str1,sql_str2];
            x = fetchmysql(sql_str,2);
        end
        %��ȡ����ָ��.ӯ������������
        %symbol���������ڣ���ƽ������ڣ�ʱ�䳤�ȣ����ͣ���ֵ
        function x = get_YingLiNengLi(var_name)
            sql_str1 = ['select symbol,publishdate,enddate,',...
                '%s from yuqerdata.yq_FdmtIndiRtnPitGet '];
            sql_str2 = [' where enddate>''2007-01-01'' ',...
                'and (secID like ''0%'' or secID like ''60%'' or secID like ''30%'')  ',... 
                'order by secID,enddate desc'];
            sql_str1 = sprintf(sql_str1,var_name);
            sql_str = [sql_str1,sql_str2];
            x = fetchmysql(sql_str,2);
        end
        %��ò���ָ�� ����FdmtDerPitGet
        %symbol���������ڣ���ƽ������ڣ�ʱ�䳤�ȣ����ͣ���ֵ
        function x = get_CaiWu_yansheng(var_name)
            sql_str1 = ['select symbol,publishdate,enddate,',...
                '%s from yuqerdata.yq_FdmtDerPitGet '];
            sql_str2 = [' where enddate>''2007-01-01'' ',...
                'and (secID like ''0%'' or secID like ''60%'' or secID like ''30%'')  ',... 
                'order by secID,enddate desc'];
            sql_str1 = sprintf(sql_str1,var_name);
            sql_str = [sql_str1,sql_str2];
            x = fetchmysql(sql_str,2);
        end
        %��ȡ����ָ��.����
        %symbol���������ڣ���ƽ������ڣ���ֵ
        function x = get_CaiWu_DanGu(var_name)
            sql_str1 = ['select symbol,publishdate,enddate,',...
                '%s from yuqerdata.yq_FdmtIndiPSPitGet '];
            sql_str2 = [' where enddate>''2007-01-01'' ',...
                'and (secID like ''0%'' or secID like ''60%'' or secID like ''30%'')  ',... 
                'order by secID,enddate desc'];
            sql_str1 = sprintf(sql_str1,var_name);
            sql_str = [sql_str1,sql_str2];
            x = fetchmysql(sql_str,2);
        end
        %��ȡ����ָ��.��Ӫ����
        %symbol���������ڣ���ƽ������ڣ���ֵ
        function x = get_CaiWu_YunYingNengLi(var_name)
            sql_str1 = ['select symbol,publishdate,enddate,',...
                '%s from yuqerdata.yq_FdmtIndiTrnovrPitGet '];
            sql_str2 = [' where enddate>''2007-01-01'' ',...
                'and (secID like ''0%'' or secID like ''60%'' or secID like ''30%'')  ',... 
                'order by secID,enddate desc'];
            sql_str1 = sprintf(sql_str1,var_name);
            sql_str = [sql_str1,sql_str2];
            x = fetchmysql(sql_str,2);
        end
        %��ȡ�����Ȳ���ָ�� �����û�з������ڣ���Ҫ�ͱ�ı����ϲ�ѯ
        %symbol���������ڣ���ƽ������ڣ���ֵ
        function x = get_CaiWu_DanJiDu(var_name)
            sql_str1 = ['select symbol,enddate,',...
                '%s from yuqerdata.yq_FdmtIndiQGet '];
            sql_str2 = [' where enddate>''2007-01-01'' ',...
                'and (secID like ''0%'' or secID like ''60%'' or secID like ''30%'')  ',... 
                'order by secID,enddate desc'];
            sql_str1 = sprintf(sql_str1,var_name);
            sql_str = [sql_str1,sql_str2];
            x = fetchmysql(sql_str,2);
        end
        %��ȡĳ��ʱ�������һ����ҵ����
        function x = get_industry_class(t)
            str_str1 = ['select ticker,industryID1 from yuqerdata.yq_industry where ',...
                'industryVersionCD=''010303'' and intodate <= ''%s'' and ',...
                '(outDate>''%s'' or outDate is null)'];
            sql_str = sprintf(str_str1,t,t);
            x = fetchmysql(sql_str,2);
        end
        function x = get_industry_class_2(t)
            str_str1 = ['select ticker,industryID1 from yuqerdata.yq_industry_sw where ',...
                'industryVersionCD=''010303'' and intodate <= ''%s'' and ',...
                '(outDate>''%s'' or outDate is null)'];
            sql_str = sprintf(str_str1,t,t);
            x = fetchmysql(sql_str,2);
        end
        function [x_st_symbol,x_st_date0,x_st_date1] = get_st_date()
            %����ST��Ϣ����
            sql_str = 'SELECT * FROM yuqerdata.st_info order by tradedate desc';
            x_st = fetchmysql(sql_str,2);
            x_st(:,1) = cellfun(@str2double,x_st(:,1),'UniformOutput',false);
            x_st_codenum = cell2mat(x_st(:,1));
            x_st_u_codenum = unique(x_st_codenum);
            x_st_data = cell(length(x_st_u_codenum),3);
            for i = 1:length(x_st_u_codenum)
                sub_x_st_data=x_st(eq(x_st_codenum,x_st_u_codenum(i)),:);
                x_st_data(i,:) = {sprintf('%0.6d',x_st_u_codenum(i)),sub_x_st_data{1,2},sub_x_st_data{end,2}};
            end
            x_st_symbol = x_st_data(:,1);
            x_st_date0 = datenum(x_st_data(:,3));
            x_st_date1 = datenum(x_st_data(:,2));
        end
        %��ȡĳ��ʱ�����ֵ
        function x = get_market_value(t)
            %��ǰ��Ʊ��ֵ
            sql_str_f1 = 'select symbol,log(marketValue) from yuqerdata.yq_dayprice where tradedate = ''%s'' and marketValue is not null and marketValue !=0 order by symbol';
            x = fetchmysql(sprintf(sql_str_f1,t),2);
        end
        %��ȡĳ��ʱ�����ͨ��ֵ
        function x = get_market_value_lt(t)
            %��ǰ��Ʊ��ֵ
            sql_str_f1 = 'select symbol,log(negMarketValue) from yuqerdata.yq_dayprice where tradedate = ''%s'' and marketValue is not null and marketValue !=0 order by symbol';
            x = fetchmysql(sprintf(sql_str_f1,t),2);
        end
        %�������
        function [x1,tref1] = filling_data(tref1,tref2,x2)
            if iscell(tref1)
                tref1 = datenum(tref1);
                tref2 = datenum(tref2);
            end
            [tref2,ia] = sort(tref2);
            x2 = x2(ia);
            T = length(tref2);
            tref1 = tref1(tref1>tref2(1));
            x1 = zeros(size(tref1));
            for i = 1:T
                if ~eq(i,T)
                    sub_ind = tref1>tref2(i) & tref1<=tref2(i+1);
                else
                    sub_ind = tref1>tref2(i);
                end
                x1(sub_ind) = x2(i);
            end
        end
        %��ʱ��������ǰ��һ��ʱ��������
        function [x1,tref1] = find_near_data(tref1,tref2,x2)
            if iscell(tref1)
                tref1 = datenum(tref1);
                tref2 = datenum(tref2);
            end
            [tref2,ia] = sort(tref2);
            x2 = x2(ia);
            T = length(tref1);
            x1 = nan(size(tref1));
            for i = 1:T
                ia = find(tref2<tref1(i),1,'last');
                if ~isempty(ia)
                    x1(i) = x2(ia);
                end
            end
            nan_ind = isnan(x1);
            x1(nan_ind) = [];
            tref1(nan_ind) = [];
        end
        function y = trans_dummy(x)
            u_x = unique(x);
            T = length(u_x);
            y = zeros(length(x),T);
            for i = 1:T
                ind = eq(x,u_x(i));
                y(ind,i) = 1;
            end
        end
        %s19 added
        %ST,*ST,PT code
        function x = get_stpt_symbol(t1)
            sql_str= 'select distinct(ticker) from yuqerdata.st_info where tradedate = ''%s''';
            x = fetchmysql(sprintf(sql_str,t1),2);
            if ~isempty(x)
                x = cellfun(@str2double,x,'UniformOutput',false);
                x = cellfun(@(x) sprintf('%0.6d',x),x,'UniformOutput',false);
            end
        end
        %��ȡ�����������Ƶ�symbol
        function x = get_time_cut_symbol(t0)
            sql_str = ['select ticker from yuqerdata.equget where listdate<''%s'' ',...
                'and listdate is not null'];
            x = fetchmysql(sprintf(sql_str,t0),2);
        end
        %��ȡĳ�յĿ���symbol
        function x = get_stop_symbol(t0)
            sql_str = ['select distinct(symbol) from yuqerdata.yq_stop_run_data where ',...
                'haltbegintime <= ''%s 09:30:00'' and haltEndTime>=''%s''' ];
            x = fetchmysql(sprintf(sql_str,t0,t0),2);
        end
		%��ȡ��Ʊ��
        function sub_symbol_pool = get_index_pool(index_pool,t_str)
            if eq(length(index_pool),6)
                sub_t = fetchmysql(sprintf(['select tradingdate from yuqerdata.IdxCloseWeightGet ',...
                    'where tradingdate < ''%s'' and ticker = ''%s'' order by tradingdate desc limit 1'],...
                                t_str,index_pool),2);
                if isempty(sub_t)
                    sub_t = fetchmysql(sprintf(['select tradingdate from yuqerdata.IdxCloseWeightGet ',...
                        'where tradingdate >= ''%s'' and ticker = ''%s''  order by tradingdate limit 1'],...
                    t_str,index_pool),2);
                end
                sub_symbol_pool = fetchmysql(sprintf(['select symbol from yuqerdata.IdxCloseWeightGet ',...
                    'where tradingdate = ''%s'' and ticker = ''%s'''],sub_t{1},index_pool),2);
            else
                if strcmpi(index_pool,'a')
                    sub_t = fetchmysql(sprintf(['select tradingdate from yuqerdata.yq_tradingdate where ',...
                        'tradingdate <=''%s'' order by tradingdate desc limit 1'],t_str),2);
                    sub_t = sub_t{1};
                    sub_symbol_pool = fetchmysql(sprintf(['select symbol from yuqerdata.yq_dayprice where ',...
                        'tradeDate = ''%s'''],sub_t),2);
                end
            end
        end
        %��ȡ�µ�����
        function tref_month = get_month_end()
            sql_str = 'select endDate from yuqerdata.yq_index_month where symbol = ''000001'' order by endDate';
            tref_month = fetchmysql(sql_str,2);
        end
        %��ȡ��ĩ����
        function tref_week = get_week_end()
            sql_str = 'select endDate from yuqerdata.yq_MktIdxwGet where ticker = ''000001'' order by endDate';
            tref_week = fetchmysql(sql_str,2);
        end
        %��ȡ����ʱ�������棬����tref*symbol����
        function [sub_r,sub_t_u,sub_symbol_u] = get_interchgPct(sub_t1,sub_t2)
            sql_str = ['select symbol,tradeDate,chgPct from yuqerdata.yq_dayprice ',...
                'where tradeDate >= ''%s'' and tradeDate <=''%s'' and chgPct is not null order by tradeDate'];
            sub_x = fetchmysql(sprintf(sql_str,sub_t1,sub_t2),2);
            %����
            sub_t_u = unique(sub_x(:,2));
            sub_symbol_u = unique(sub_x(:,1));
            T_sub_t_u = length(sub_t_u);
            T_sub_symbol = length(sub_symbol_u);
            sub_r = zeros(T_sub_symbol,T_sub_t_u);
            for j = 1:T_sub_t_u
                sub_sub_x = sub_x(strcmp(sub_x(:,2),sub_t_u(j)),[1,3]);
                [~,ia,ib] = intersect(sub_symbol_u,sub_sub_x(:,1));
                sub_r(ia,j) = cell2mat(sub_sub_x(ib,2));
            end
        end
        %��ȡ��ʼʱ��
        function t0 = get_table_end_date(tn,var_name)
            if nargin<2
                var_name = 'tradingdate';
            end
            sql_str = 'select %s from %s order by %s desc limit 1';
            t0 = fetchmysql(sprintf(sql_str,var_name,tn,var_name),2);
            if ~isempty(t0)
                t0= t0{1};
            end
        end
    end
    
    
end