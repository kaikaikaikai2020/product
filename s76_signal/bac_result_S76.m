classdef bac_result_S76 < handle
    
    methods
        function get_all_results(obj)
            %升级股票池、信号
            dos('python bac_toolS76.py')
            dos('python bac_toolS76P2.py')
            %sprintf('开始更新S75P1HK数据')
            %升级曲线
            obj.bac_S76P1()
            obj.get_group_symbol_S76P1()
            %计算结果并输出
            [H1,re1,info1] = obj.bac_fig_76();
            [H2,re2,info2] = obj.bac_fig_76P2();
            [H3,re3,info3] = obj.bac_re_csi();
            [H4,re4,info4] = bac_re_iso();
%             H2 = [];
%             re2 = [];
%             info2 = [];
            key_str = 'S76综合项目';
            file_name = sprintf('%s%s',key_str,datestr(now,'yyyy-mm-dd'));
            pn0 = fullfile(pwd,'计算结果');
            if ~exist(pn0,'dir')
                mkdir(pn0);
            end
            obj_wd = wordcom(fullfile(pn0,sprintf('%s.doc',file_name)));
            H = [H1,H2,H3,H4];
            for i = 1:length(H)
                obj_wd.pasteFigure(H(i),' ');
            end
            obj_wd.CloseWord();
            yc = [re1,re2,re3,re4];
            yt = [info1,info2,info3,info4];
            %yt = [yt{:}];
            sta_re = obj.curve_static_batch(yc,yt);
            xlstocsv_adair(fullfile(pn0,sprintf('%s.xlsx',file_name)),sta_re) 
        end
    end
    
    methods(Static)
        
        %get_group_symbol_S73P1
        function get_group_symbol_S76P1()
            key_str = 'S76均衡估值因子';    
            
            symbol_pool0 = {'A','000001','000300','000905'};
            
            m_pool = {'BET_c', 'BET_c_STD', 'PEG_c', 'PEG_c_STD'};
            tN = 's37.symbol_pool_s76p1';
            
            tN_symbol =  'zx02.symbol_pool_s76';
            var_info = {'taskID','mID','tradeDate','poolID','pool_l','pool_s'};

            re_f0 = cell(size(symbol_pool0));
            for pool_Id = 1:length(symbol_pool0)
                symbol_pool = symbol_pool0{pool_Id};
                re= cell(size(m_pool));
                for i0 = 1:length(m_pool)        
                    m_id = m_pool{i0};
                    t0_start = 'select tradeDate from %s where poolID="%s" and mID="%s" order by tradeDate desc limit 1';
                    t0_start = fetchmysql(sprintf(t0_start,tN_symbol,symbol_pool,m_id),2);
                    t0_ini = '2003-01-01';
                    if isempty(t0_start)
                        t0_start = t0_ini;
                    else
                        t0_start = t0_start{1};
                    end
                    sql_str = 'select distinct(tradeDate) from %s where tradeDate>"%s" and %s is not null order by tradeDate';
                    p_t = fetchmysql(sprintf(sql_str,tN,t0_start,m_id),2);
                    %间隔时间
                    id_t0=1;                    
                    T = length(p_t);
                    sql_str0 = 'select ticker,%s from %s where tradeDate = "%s" and %s is not null';
                    r = cell(T,1);
                    parfor i = id_t0:T
                        sub_t1 = p_t{i};                        
                        sub_ticker =fetchmysql(sprintf(sql_str0,m_id,tN,sub_t1,m_id),2);        
                        %股票池限制
                        if eq(length(symbol_pool),6)
                            sub_ticker_limit = yq_methods.get_index_pool(symbol_pool,sub_t1);
                        else
                            sub_ticker_limit = sub_ticker(:,1);
                        end
                        [~,ia] = intersect(sub_ticker(:,1),sub_ticker_limit);
                        sub_ticker = sub_ticker(ia,:);
                        f = cell2mat(sub_ticker(:,2));
                        [~,ia]=sort(f,'descend');
                        sub_ticker = sub_ticker(ia,:);                        
                        %拆分;
                        ticker_num = floor(length(f)/10);
                        sub_ticker_l = sub_ticker(1:ticker_num,1);
                        sub_ticker_s = sub_ticker(end-ticker_num+1:end,1);
                        sub_ticker_l = strjoin(sub_ticker_l,',');
                        sub_ticker_s = strjoin(sub_ticker_s,',');
                        r{i}={'S76P1',m_id,sub_t1,symbol_pool,sub_ticker_l,sub_ticker_s}';
                        sprintf('%s\n%s-%s %d-%d',key_str,symbol_pool,m_id,i,T)
                    end
                    re{i0} = [r{:}];
                end
                re = [re{:}];
                re_f0{pool_Id} = re;
            end
            re_f0=[re_f0{:}]';
            if ~isempty(re_f0)
                datainsert_adair(tN_symbol,var_info,re_f0)
            end
        end
        
        function bac_S76P1()
            key_str = 'S76均衡估值因子';            
            symbol_pool0 = {'A','000001','000300','000905'};
            %symbol_pool0 = {'A'};
            m_pool = {'BET_c', 'BET_c_STD', 'PEG_c', 'PEG_c_STD'};
            tN = 's37.symbol_pool_s76p1';

            tN_chg = 's37.s76_return';
            var_info = {'pID','mID','tradeDate','long_r','short_r','ls_r'};

            re_f0 = cell(size(symbol_pool0));
            for pool_Id = 1:length(symbol_pool0)
                symbol_pool = symbol_pool0{pool_Id};
                re= cell(size(m_pool));
                for i0 = 1:length(m_pool)        
                    m_id = m_pool{i0};
                    t0_start = 'select tradeDate from %s where pID="%s" and mID="%s" order by tradeDate desc limit 1';
                    t0_start = fetchmysql(sprintf(t0_start,tN_chg,symbol_pool,m_id),2);
                    t0_ini = '2009-11-01';
                    if isempty(t0_start)
                        t0_start = t0_ini;
                    else
                        t0_start = t0_start{1};
                    end
                    sql_str = 'select distinct(tradeDate) from %s where "%s" is not null order by tradeDate';
                    p_t = fetchmysql(sprintf(sql_str,tN,m_id),2);
                    fee = 3/1000;
                    %fee = 0;
                    %间隔时间
                    p_t_id = datenum(p_t);
                    id_t0=find(p_t_id<=datenum(t0_start),1,'last');
                    tref = yq_methods.get_tradingdate(t0_start);
                    tref = tref(2:end);
                    if isempty(tref)
                        sprintf('%s %s-%s已经是最新:%s',key_str,symbol_pool,m_id,t0_start)
                        continue
                    end
                    %tref_num = datenum(tref);
                    T = length(p_t);
                    %sql_str = ['select symbol,tradeDate,chgPct from yuqerdata.yq_dayprice where symbol ',...
                    %    'in (%s) and tradeDate > "%s" and tradeDate<= "%s" order by tradeDate'];
                    sql_str = ['select symbol,tradeDate,chgPct from yuqerdata.yq_dayprice where ',...
                        'tradeDate > "%s" and tradeDate<= "%s" order by tradeDate'];
                    
                    sql_str0 = 'select ticker,%s from %s where tradeDate = "%s" and %s is not null';
                    r = cell(T,1);
                    tref2 = cell(T,1);
                    parfor i = id_t0:T
                        sub_t1 = p_t{i};
                        if i < T
                            sub_t2 = p_t{i+1};
                        else
                            sub_t2 = tref{end};
                        end
                        %sub_ind = strcmp(tref_pub,p_t(i));
                        sub_ticker =fetchmysql(sprintf(sql_str0,m_id,tN,sub_t1,m_id),2);        
                        %股票池限制
                        if eq(length(symbol_pool),6)
                            sub_ticker_limit = yq_methods.get_index_pool(symbol_pool,sub_t1);
                        else
                            sub_ticker_limit = sub_ticker(:,1);
                        end
                        [~,ia] = intersect(sub_ticker(:,1),sub_ticker_limit);
                        sub_ticker = sub_ticker(ia,:);
                        f = cell2mat(sub_ticker(:,2));
                        [~,ia]=sort(f,'descend');
                        sub_ticker = sub_ticker(ia,:);
                        %sub_ticker1 = sub_ticker(ia(1:floor(length(sub_ticker)*0.1)),1);
                        %sub_ticker1 = sub_ticker(ia(1:ticker_num),1);
                        %sub_ticker2 = sub_ticker(ia(end-ticker_num+1:end),1);
                        %sub_ticker = [sub_ticker1;sub_ticker2];

                        %sub_info = sprintf('"%s"',strjoin(sub_ticker,'","'));    
                        sub_r = fetchmysql(sprintf(sql_str,sub_t1,sub_t2),3);
                        if isempty(sub_r)
                            continue
                        end
                        sub_tref = unique(sub_r.tradeDate);
                        tref2{i} = sub_tref';
                        sub_r = unstack(sub_r,'chgPct','tradeDate');
                        sub_r = table2cell(sub_r)';
                        sub_symbol = sub_r(1,:);
                        
                        [~,ia,ib] = intersect(sub_ticker(:,1),sub_symbol,'stable');
                        
                        sub_r0 = cell(size(sub_r,1),length(sub_ticker(:,1)));
                        sub_r0(:,ia) = sub_r(:,ib);
                        
                        %sub_r0 = sub_r(:,ib);
                        sub_r = cell2mat(sub_r0(2:end,:));
                        sub_r(isnan(sub_r)) = 0;
                        %拆分;
                        ticker_num = floor(length(sub_r(1,:))/10);
                        sub_r1 = sub_r(:,1:ticker_num);
                        sub_r2 = sub_r(:,end-ticker_num+1:end);
                        sub_r1 = mean(sub_r1,2);
                        sub_r2 = mean(sub_r2,2);
                        sub_r3 = sub_r1-sub_r2;
                        
                        sub_r_re = [sub_r1,sub_r2,sub_r3];
                        sub_r_re(1,:) = sub_r_re(1,:) -fee/2;
                        sub_r_re(end,:) = sub_r_re(end,:) -fee/2;
                        r{i}=sub_r_re';
                        sprintf('%s\n%s-%s %d-%d',key_str,symbol_pool,m_id,i,T)
                    end
                    tref2 = [tref2{:}]';
                    y = [r{:}]';
                    sub_re = [tref2,tref2,tref2,num2cell(y)];
                    sub_re(:,1) = {symbol_pool};
                    sub_re(:,2) = {m_id};
                    [~,ia] = intersect(sub_re(:,3),tref);
                    re{i0} = sub_re(ia,:)';
                end
                re = [re{:}];
                re_f0{pool_Id} = re;
            end
            re_f0=[re_f0{:}]';
            if ~isempty(re_f0)
                datainsert_adair(tN_chg,var_info,re_f0)
            end
        end
        
                
        function [H,Y,TN] =bac_fig_76()
            method_key = {'BET_c', 'BET_c_STD', 'PEG_c', 'PEG_c_STD'};
            method_info = containers.Map(method_key,...
                {'静态BET','动态BET','静态PEG','动态PEG'});
            tN = 's37.s76_return';
            index_symbol = {'A','000001','000905','000300'};
            index_info = containers.Map(index_symbol,{'全A股票池','上证综指股票池','500股票池','300股票池'});
            %index_symbol = {'A'};
            %index_info = containers.Map(index_symbol,{'全A股票池'});
            sql_str1 = 'select * from %s where mID = "%s" and pID = "%s" order by tradeDate';
            H=[];
            Y={};
            TN ={};
            for m_num = 1:length(method_key)

                m_id = method_key;
                m_id = m_id{m_num};
                for j = 1:length(index_symbol)
                    p_id=index_symbol{j};
                    x = fetchmysql(sprintf(sql_str1,tN,m_id,p_id),2);
                    %tref = x(:,3);
                    %[tref,ia,ib] = intersect(tref,index_chg{j}(:,1));
                    %x = cell2mat(x(ia,end-2:end));
                    %y = cell2mat(index_chg{j}(ib,2));
                    tref = x(:,3);
                    x = cell2mat(x(:,end-2:end));
                    x(1,:) = 0;
                    %y(1) = 0;

                    %yc1 = cumprod(1+x(:,1));
                    %yc2 = cumprod(1+x-y);        
                    %yc2 = cumprod(1+x(:,2));
                    yc3 = cumprod(1+x(:,3));
                    t_str = sprintf('%s-%s',method_info(m_id),index_info(p_id));        
                    %t_str1 = sprintf('%s-long',t_str);
                    %t_str2 = sprintf('%s-short',t_str);        
                    t_str3 = sprintf('%s-long-short',t_str);
                    t_str3 = replace(t_str3,'_','-');

                    %h1 = figure_S53(yc1,tref,t_str1);
                    %h2 = figure_S53(yc2,tref,t_str2);
                    h3 = figure_S53(yc3,tref,t_str3);

                    H = cat(2,H,[h3]);
                    Y = cat(2,Y,{yc3});
                    TN = cat(2,TN,{t_str3});
                end
            end
        end
        %指数成分股选股结果
        function [H,Y,info] = bac_re_iso()
            key_str = sprintf('S76 趋势策略成分股组合%s',datestr(now,'yyyymmdd'));
            %forex-day,as51,topix
            index_pool = {'as51', 'topix', 'twse', 'kosdaq', ...
                'kospi', 'msci', 'ndx', 'nifty', 'nky', 'RTY', 'set50', 'sx5e', 'ukx', 'xin9i'};
            
            
            index_info = containers.Map({'as51','topix','twse','hsce','hkggt','kosdaq',...
                'kospi', 'msci', 'ndx', 'nifty', 'nky', 'RTY', 'set50', 'sx5e',...
                       'ukx', 'xin9i'},{'AS51','TPX','TWSE','HSCEI','hkggt','KOSDAQ',...
                       'KOSPI2','TAMSCI','NDX','NIFTY','NKY','RTY','SET50','SX5E',...
                      'UKX','XIN9I'});
            T_index_pool = length(index_pool);
            Y = cell(1,T_index_pool);
            H = zeros(1,T_index_pool);
            info = Y;
            for i = 1:T_index_pool
                index0 =index_pool{i};
                sql_str = sprintf('select ticker,tradeDate,CHGPct,sig  from S37.s76_%s order by tradeDate',index0);
                x0 = fetchmysql(sql_str,2);
                x = x0(:,[1,2,4]);
                x = cell2table(x,'VariableNames',{'s','t','v'});
                ind = unstack(x,'v','s');
                ind.Properties.RowNames = ind.t;
                ind.t = [];
                pn0 = write_signal_check();
                file_name = sprintf('signal_%s_%s.csv',index0,x0{end,2});
                writetable(ind,fullfile(pn0,file_name),'WriteRowNames',true)
                x = x0(:,[1,2,3]);
                x = cell2table(x,'VariableNames',{'s','t','v'});
                X = unstack(x,'v','t');

                tref = X.Properties.VariableNames;
                tref = cellfun(@(x) x(2:end),tref,'UniformOutput',false);
                tref = tref(2:end);
                X = table2cell(X); 
                symbols = X(:,1);
                X = cell2mat(X(:,2:end))';
                X(isnan(X)) = 0;

                ind0 = eq(sum(abs(X)),0);
                X(:,ind0) = [];
                symbols(ind0) = [];

                t_str = cellfun(@(x) [x(1:4),x(6:7),x(9:10)],tref,'UniformOutput',false)';
                T = length(t_str);

                title_str = index0;
                title_str(strfind(title_str,'_')) = '-';

                sub_symbols = symbols;

                [~,ia] = intersect(symbols,sub_symbols);
                sub_Y = X(:,ia);                
                sub_r = mean(sub_Y,2);
                
                sql_tmp = 'select tradeDate,closePrice from data_pro.main_index_s68 where index_id = "%s" and ticker = "%s" order by tradeDate';
                sub_r2 = fetchmysql(sprintf(sql_tmp,index0,index_info(index0)),2);
                
                sub_tref2 = cellfun(@(x) replace(x,'-',''),sub_r2(:,1),'UniformOutput',false);
                sub_r2 = cell2mat(sub_r2(:,2));
                sub_r2(2:end) = sub_r2(2:end)./sub_r2(1:end-1)-1;
                sub_r2(1) = 0;                
                [t_str,ia,ib] = intersect(t_str,sub_tref2);
                sub_r = sub_r(ia);
                sub_r2 = sub_r2(ib);
                
                r_c = cumprod(1+[sub_r,sub_r2,sub_r-sub_r2]);
                leg_str = {'S76P1',index0,['S76P1-',index0]};
                h = bacFigure(r_c,t_str,title_str,leg_str,1);
                
                Y{i} = {r_c(:,1),r_c(:,2),r_c(:,3)};
                H(i) = h;
                info{i} = leg_str;
                sprintf('%s-%d-%d',key_str,i,T_index_pool)
            end
            
            Y = [Y{:}];
            info = [info{:}];
            
        end
        function [H,Y,TN] =bac_fig_76P2()
            %后续修改这个表格的名称就可以了
            tn = 'parapool.s76p2_signal_signal1';
            info = sprintf('select distinct(concat(dtype,"-",ticker)) from %s where ticker not in ("SM","RS","PM","TS","LR","USDCNH")',tn);
            info = fetchmysql(info,2);
            info = cellfun(@(x) strsplit(x,'-'),info,'UniformOutput',false);
            info1 = cellfun(@(x) x{1},info,'UniformOutput',false);
            info2 = cellfun(@(x) x{2},info,'UniformOutput',false);
            
            sql_str1 = 'select tradeDate,CHGPct,sig from %s where dtype = "%s" and ticker = "%s" order by tradeDate';
            H=[];
            Y={};
            TN ={};
            for m_num = 1:length(info1)
                x = fetchmysql(sprintf(sql_str1,tn,info1{m_num},info2{m_num}),2);
                tref = x(:,1);
                y1 = cell2mat(x(:,2));
                y2 = cell2mat(x(:,3));
                
                y2(2:end) = y2(1:end-1);
                y2(1) = 0;
                
                x = y1.*y2;
                yc3 = [cumprod(1+y1),cumprod(1+x)];
                if yc3(end,1)>yc3(end,2)
                    continue
                end
                t_str = sprintf('%s-%s',info1{m_num},info2{m_num});        
                t_str3 = replace(t_str,'_','-');
                leg_str = {t_str3,[t_str3,'-择时']};
                h3 = figure_S53(yc3,tref,t_str3,leg_str);

                H = cat(2,H,h3);
                Y = cat(2,Y,{yc3(:,1),yc3(:,2)});
                TN = cat(2,TN,leg_str);
            end
        end
        function [H,Y,info] = bac_re_csi()
            pn0 = write_signal_check();
            key_str = sprintf('S76csi补充%s',datestr(now,'yyyymmdd'));            
            sql_str = 'select ticker,tradeDate,CHGPct*sig1-fee,sig  from parapool.s76_csi_stock order by tradeDate';
            x0 = fetchmysql(sql_str,2);
            x = x0(:,[1,2,4]);
            x = cell2table(x,'VariableNames',{'s','t','v'});
            pool_w = unstack(x,'v','s');
            pool_w.Properties.RowNames = pool_w.t;
            pool_w.t = [];                        
            %file_name = sprintf('signal_%s_%s.csv','csi',x0{end,2});
            %writetable(ind,fullfile(pn0,file_name),'WriteRowNames',true)
            %symbols = unique(x(:,1));
            %tref = unique(x(:,2));
            x = x0(:,[1,2,3]);
            x = cell2table(x,'VariableNames',{'s','t','v'});
            X = unstack(x,'v','t');

            tref = X.Properties.VariableNames;
            tref = cellfun(@(x) x(2:end),tref,'UniformOutput',false);
            tref = tref(2:end);
            X = table2cell(X); 
            symbols = X(:,1);
            X = cell2mat(X(:,2:end))';
            X(isnan(X)) = 0;

            ind0 = eq(sum(abs(X)),0);
            X(:,ind0) = [];
            symbols(ind0) = [];

            t_str = cellfun(@(x) [x(1:4),x(6:7),x(9:10)],tref,'UniformOutput',false);
            T = length(t_str);

            symbol_pool_all = {   'a',    '000905','000300','000852'};
            symbol_pool_info = {'全市场','中证500','沪深300','中证1000'};
            sql_index = 'select tradeDate,CHGPct from yuqerdata.yq_index where symbol = "%s" order by tradeDate';
            
            
            T_index_pool = length(symbol_pool_all);
            Y = cell(1,T_index_pool);
            H = zeros(1,T_index_pool);
            info = Y;
            symbols_comp = cell(T_index_pool,1);
            for i = 1:T_index_pool
                sub_index = symbol_pool_all{i};
                title_str = symbol_pool_info{i};
                title_str(strfind(title_str,'_')) = '-';
                sub_symbols = yq_methods.get_index_pool(sub_index,datestr(now,'yyyy-mm-dd'));
                tmp = cellfun(@(x) ['x',x],sub_symbols,'UniformOutput',false);
                [~,ia] = intersect(pool_w.Properties.VariableNames,tmp);
                sub_pool_w = pool_w(:,ia);
                file_name = sprintf('signal_%s_csi_%s.csv',title_str,x0{end,2});
                writetable(sub_pool_w,fullfile(pn0,file_name),'WriteRowNames',true)
            
                symbols_comp{i} = sub_symbols;
                [~,ia] = intersect(symbols,sub_symbols);
                sub_Y = X(:,ia);            
                sub_r = mean(sub_Y,2);
                
                if length(sub_index)<6
                    tmp = '000001';
                else
                    tmp = sub_index;
                end
                sub_r2 = sprintf(sql_index,tmp);
                sub_r2 = fetchmysql(sub_r2,2);
                sub_r2(:,1) = cellfun(@(x) replace(x,'-',''),sub_r2(:,1),'UniformOutput',false);
                [t_str,ia,ib] = intersect(t_str,sub_r2(:,1));
                sub_r = [sub_r(ia),cell2mat(sub_r2(ib,2)),sub_r(ia)-cell2mat(sub_r2(ib,2))];
                r_c = cumprod(1+sub_r);
                
                leg_str = {'S76P1',sub_index,['S76P1-',sub_index]};
                h = bacFigure(r_c,t_str,title_str,leg_str,1);
                
                Y{i} = {r_c(:,1),r_c(:,2),r_c(:,3)};
                H(i) = h;
                info{i} = leg_str;
                
            end
            Y = [Y{:}];
            info = [info{:}];
        end
        function sta_re = curve_static_batch(yc,title_str)
            if ~iscell(yc)
                temp = cell(size(yc(1,:)));
                for i = 1:size(yc,2)
                    temp{i} = yc(:,i);
                end
                yc = temp;
            end
            sta_re = cell(size(yc));
            for i = 1:length(sta_re)
                [v0,v_str0] = curve_static(yc{i},[]);
                [v,v_str] = ad_trans_sta_info(v0,v_str0);
                if eq(i,1)
                    sub_re = [[{''};v_str'],[title_str(i);v']];
                else
                    sub_re = [title_str(i);v'];
                end
                sta_re{i} = sub_re;
            end
            sta_re = [sta_re{:}]';
        end
    end
    
end

function h = figure_S53(y_re,tref,t_str1,leg_str)
if nargin <4
    leg_str = [];
end
t_str = cellfun(@(x) strjoin(strsplit(x,'-'),''),tref,'UniformOutput',false);
T = length(t_str);
h = figure;
plot(y_re,'LineWidth',2);
set(gca,'xlim',[0,T]);
set(gca,'XTick',floor(linspace(1,T,15)));
set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
set(gca,'XTickLabelRotation',90)    
setpixelposition(h,[223,365,1345,420]);
if ~isempty(leg_str)
legend(leg_str,'Location','best')
end
box off
title(t_str1)
end
function pn0 = write_signal_check()
pn0 = fullfile(pwd,'para_pool_S76');
if ~exist(pn0,'dir')
    mkdir(pn0);
end
end