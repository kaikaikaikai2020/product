classdef bac_result_S53 < handle
    
    methods
        function get_all_results(obj)
            %升级股票池
            obj.update_signal_p1();
            obj.update_signal_p2();
            %升级曲线
            obj.bac_S53P1();
            obj.bac_S53P2();
            %计算结果并输出
            [H1,re1,info1] = obj.bac_fig_53P1();
            [H2,re2,info2] = obj.bac_fig_53P2();
            
            key_str = 'S51 基金重仓股独门股等策略';
            file_name = sprintf('%s%s',key_str,datestr(now,'yyyy-mm-dd'));
            pn0 = fullfile(pwd,'计算结果');
            if ~exist(pn0,'dir')
                mkdir(pn0);
            end
            obj_wd = wordcom(fullfile(pn0,sprintf('%s.doc',file_name)));
            H = [H1,H2];
            for i = 1:length(H)
                obj_wd.pasteFigure(H(i),' ');
            end
            obj_wd.CloseWord();
            yc = [re1,re2];
            yt = [info1,info2];
            %yt = [yt{:}];
            sta_re = obj.curve_static_batch(yc,yt);
            xlstocsv_adair(fullfile(pn0,sprintf('%s.xlsx',file_name)),sta_re) 
        end
    end
    
    methods(Static)
        function update_signal_p1()
            sprintf('开始更新S53P1(MHKQ因子择时模型在A股中的运用)因子数据，耗时10分钟左右')
            dos('python S53_bac_tool1.py')
            sprintf('完成更新S53P1(MHKQ因子择时模型在A股中的运用)因子数据，耗时10分钟左右')
        end
        function update_signal_p2()
            sprintf('开始更新S53P2(细分行业下的多因子选股模型)因子数据，耗时10分钟左右')
            dos('python S53_bac_tool2_update2.py')
            sprintf('完成更新S53P2(细分行业下的多因子选股模型)因子数据，耗时10分钟左右')
        end        
        function bac_S53P1()
            key_str = 'S53P1 MHKQ因子择时模型在A股中的运用 收益计算';
            symbol_pool0 = {'A','000300','000905'};
            m_pool = {'factor_MHKQ','factor_IC'};
            tN = 'symbol_pool_s53P1';

            tN_chg = 'S37.S53_return';
            var_info = {'pID','mID','tradeDate','r'};

            re_f0 = cell(size(symbol_pool0));
            for pool_Id = 1:length(symbol_pool0)
                symbol_pool = symbol_pool0{pool_Id};
                re= cell(size(m_pool));
                for i0 = 1:length(m_pool)        
                    m_id = m_pool{i0};
                    t0_start = 'select max(tradeDate) from %s where pID="%s" and mID="%s"';
                    t0_start = fetchmysql(sprintf(t0_start,tN_chg,symbol_pool,m_id),2);
                    if isempty(t0_start)
                        t0_start = '2011-02-01';
                    else
                        t0_start = t0_start{1};
                    end
                    sql_str = 'select distinct(tradeDate) from S37.%s where tradeDate >="2011-02-01" order by tradeDate';
                    p_t = fetchmysql(sprintf(sql_str,tN),2);
                    fee = 2/1000;
                    %fee = 0;
                    %间隔时间
                    p_t_id = datenum(p_t);
                    id_t0=find(p_t_id<=datenum(t0_start),1,'last');
                    %tref = yq_methods.get_tradingdate(p_t{id_t0});
                    tref = yq_methods.get_tradingdate(t0_start);
                    tref = tref(2:end);
                    if isempty(tref)
                        sprintf('%s %s-%s已经是最新:%s',key_str,symbol_pool,m_id,t0_start)
                        continue
                    end
                    tref_num = datenum(tref);
                    T = length(p_t);
                    sql_str = ['select symbol,tradeDate,chgPct from yuqerdata.yq_dayprice where symbol ',...
                        'in (%s) and tradeDate > "%s" and tradeDate<= "%s" order by tradeDate'];
                    sql_str0 = 'select code,%s from S37.%s where tradeDate = "%s"';
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
                        sub_ticker =fetchmysql(sprintf(sql_str0,m_id,tN,sub_t1),2);        
                        %股票池限制
                        if eq(length(symbol_pool),6)
                            sub_ticker_limit = yq_methods.get_index_pool(symbol_pool,sub_t1);
                        else
                            sub_ticker_limit = sub_ticker(:,1);
                        end
                        [~,ia] = intersect(sub_ticker(:,1),sub_ticker_limit);
                        sub_ticker = sub_ticker(ia,:);
                %         %st pt delete
                %         symbol_stpt = yq_methods.get_stpt_symbol(sub_t1);
                %         [~,ia] = intersect(sub_ticker(:,1),symbol_stpt);
                %         sub_ticker(ia,:) = [];
                %         %上市限制
                %         symbol_old = list_datelimit(sub_t1,250);
                %         [~,ia] = intersect(sub_ticker(:,1),symbol_old);
                %         sub_ticker = sub_ticker(ia,:);
                        f = cell2mat(sub_ticker(:,2));
                        [~,ia]=sort(f,'descend');

                        %sub_ticker1 = sub_ticker(ia(1:floor(length(sub_ticker)*0.1)),1);
                        sub_ticker1 = sub_ticker(ia(1:100),1);

                        sub_ticker2 = sub_ticker(ia(end-floor(length(sub_ticker)*0.1)+1:end),1);
                        sub_ticker = [sub_ticker1;sub_ticker2];

                        sub_info = sprintf('"%s"',strjoin(sub_ticker,'","'));    
                        sub_r = fetchmysql(sprintf(sql_str,sub_info,sub_t1,sub_t2),3);
                        if isempty(sub_r)
                            continue
                        end
                        sub_tref = unique(sub_r.tradeDate);
                        tref2{i} = sub_tref';
                        sub_r = unstack(sub_r,'chgPct','tradeDate');
                        sub_r = table2cell(sub_r)';
                        sub_symbol = sub_r(1,:);

                        [sub_symbol,ia,ib] = intersect(sub_symbol,sub_ticker,'stable');
                        sub_r0 = cell(size(sub_r,1),length(sub_ticker));
                        sub_r0(:,ib) = sub_r(:,ia);

                        sub_r = cell2mat(sub_r0(2:end,:));
                        sub_r(isnan(sub_r)) = 0;
                        %拆分
                        ia = length(sub_ticker1);
                        sub_r1 = sub_r(:,1:ia);
                        sub_r2 = sub_r(:,ia+1:end);
                        sub_r1 = mean(sub_r1,2);
                        sub_r2 = mean(sub_r2,2);
                        sub_r_re = sub_r1-0;
                        sub_r_re(end) = sub_r_re(end)-fee;
                        %sub_r(end) = sub_r(1)-fee;
                        %sub_r(1) = -fee;
                        %sub_r_re = mean(sub_r,2);
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
        function bac_S53P2()
            key_str = 'S53P2 细分行业下的多因子选股模型 收益计算';
            symbol_pool0 = {'A','000300','000905'};
            m_pool = {'S53P2basic','S53P2final'};
            tN = 'symbol_pool_s53P2';
            tN_chg = 'S37.S53_return';
            var_info = {'pID','mID','tradeDate','r'};

            re_f0 = cell(size(symbol_pool0));
            for pool_Id = 1:length(symbol_pool0)
                symbol_pool = symbol_pool0{pool_Id};
                re= cell(size(m_pool));
                for i0 = 1:length(m_pool)        
                    m_id = m_pool{i0};
                    t0_start = 'select tradeDate from %s where pID="%s" and mID="%s" order by tradeDate desc limit 1';
                    t0_start = fetchmysql(sprintf(t0_start,tN_chg,symbol_pool,m_id),2);
                    if isempty(t0_start)
                        t0_start = '2011-01-01';
                    else
                        t0_start = t0_start{1};
                    end
                    sql_str = 'select distinct(tradeDate) from S37.%s  order by tradeDate';
                    p_t = fetchmysql(sprintf(sql_str,tN),2);
                    fee = 2/1000;
                    %fee = 0;
                    %间隔时间
                    p_t_id = datenum(p_t);
                    id_t0=find(p_t_id<=datenum(t0_start),1,'last');
                    %tref = yq_methods.get_tradingdate(p_t{id_t0});
                    tref = yq_methods.get_tradingdate(t0_start);
                    tref = tref(2:end);
                    if isempty(tref)
                        sprintf('%s %s-%s已经是最新:%s',key_str,symbol_pool,m_id,t0_start)
                        continue
                    end
                    tref_num = datenum(tref);
                    T = length(p_t);
                    sql_str = ['select symbol,tradeDate,chgPct from yuqerdata.yq_dayprice where symbol ',...
                        'in (%s) and tradeDate > "%s" and tradeDate<= "%s" order by tradeDate'];
                    sql_str0 = 'select ticker,f_val from S37.%s where tradeDate = "%s" and mID="%s"';
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
                        sub_ticker =fetchmysql(sprintf(sql_str0,tN,sub_t1,m_id),2);        
                        %股票池限制
                        if eq(length(symbol_pool),6)
                            sub_ticker_limit = yq_methods.get_index_pool(symbol_pool,sub_t1);
                        else
                            sub_ticker_limit = sub_ticker(:,1);
                        end
                        [~,ia] = intersect(sub_ticker(:,1),sub_ticker_limit);
                        sub_ticker = sub_ticker(ia,:);
                %         %st pt delete
                %         symbol_stpt = yq_methods.get_stpt_symbol(sub_t1);
                %         [~,ia] = intersect(sub_ticker(:,1),symbol_stpt);
                %         sub_ticker(ia,:) = [];
                %         %上市限制
                %         symbol_old = list_datelimit(sub_t1,250);
                %         [~,ia] = intersect(sub_ticker(:,1),symbol_old);
                %         sub_ticker = sub_ticker(ia,:);
                        f = cell2mat(sub_ticker(:,2));
                        [~,ia]=sort(f,'descend');

                        %sub_ticker1 = sub_ticker(ia(1:floor(length(sub_ticker)*0.1)),1);
                        sub_ticker1 = sub_ticker(ia(1:100),1);

                        sub_ticker2 = sub_ticker(ia(end-floor(length(sub_ticker)*0.1)+1:end),1);
                        sub_ticker = [sub_ticker1;sub_ticker2];

                        sub_info = sprintf('"%s"',strjoin(sub_ticker,'","'));    
                        sub_r = fetchmysql(sprintf(sql_str,sub_info,sub_t1,sub_t2),3);
                        if isempty(sub_r)
                            continue
                        end
                        sub_tref = unique(sub_r.tradeDate);
                        tref2{i} = sub_tref';
                        sub_r = unstack(sub_r,'chgPct','tradeDate');
                        sub_r = table2cell(sub_r)';
                        sub_symbol = sub_r(1,:);

                        [sub_symbol,ia,ib] = intersect(sub_symbol,sub_ticker,'stable');
                        sub_r0 = cell(size(sub_r,1),length(sub_ticker));
                        sub_r0(:,ib) = sub_r(:,ia);

                        sub_r = cell2mat(sub_r0(2:end,:));
                        sub_r(isnan(sub_r)) = 0;
                        %拆分
                        ia = length(sub_ticker1);
                        sub_r1 = sub_r(:,1:ia);
                        sub_r2 = sub_r(:,ia+1:end);
                        sub_r1 = mean(sub_r1,2);
                        sub_r2 = mean(sub_r2,2);
                        sub_r_re = sub_r1-0;
                        sub_r_re(end) = sub_r_re(end)-fee;
                        %sub_r(end) = sub_r(1)-fee;
                        %sub_r(1) = -fee;
                        %sub_r_re = mean(sub_r,2);
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
        
        function [H,Y,TN] =bac_fig_53P1()
            method_key = {'factor_IC','factor_MHKQ','S53P2basic','S53P2final'};
            method_info = containers.Map(method_key,...
                {'S53P1IC加权','S53P1MHKQ加权','S53P2全行业IC加权','S53P2分行业IC加权'});
            tN = 'S37.S53_return';
            %var_info = {'method_info','tradingdate','r','r300','r500','dr300','dr500'};
            index_symbol = {'A','000905','000300'};
            index_info = containers.Map(index_symbol,{'A股股票池','500股票池','300股票池'});
            %index result
            sql_str = 'select tradeDate,CHGPct from yuqerdata.yq_index where symbol = "%s" order by tradeDate';
            index_chg = cell(size(index_symbol));
            for i = 1:length(index_symbol)
                if strcmp(index_symbol{i},'A')
                    temp = '000985';
                else
                    temp = index_symbol{i};
                end
                index_chg{i} = fetchmysql(sprintf(sql_str,temp),2);
            end

            sql_str1 = 'select * from %s where mID = "%s" and pID = "%s" order by tradeDate';

            H=[];
            Y={};
            TN ={};
            for m_num = 3:length(method_key)

                m_id = method_key;
                m_id = m_id{m_num};
                for j = 2:length(index_symbol)
                    p_id=index_symbol{j};
                    x = fetchmysql(sprintf(sql_str1,tN,m_id,p_id),2);
                    tref = x(:,3);
                    [tref,ia,ib] = intersect(tref,index_chg{j}(:,1));
                    x = cell2mat(x(ia,end));
                    y = cell2mat(index_chg{j}(ib,2));
                    x(1) = 0;
                    y(1) = 0;

                    yc1 = cumprod(1+x);
                    yc2 = cumprod(1+x-y);        
                    t_str = sprintf('%s-%s',method_info(m_id),index_info(p_id));        
                    t_str1 = sprintf('%s-做多',t_str);
                    t_str2 = sprintf('%s-指数对冲',t_str);        

                    h1 = figure_S53(yc1,tref,t_str1);
                    h2 = figure_S53(yc2,tref,t_str2);

                    H = cat(2,H,[h1,h2]);
                    Y = cat(2,Y,{yc1,yc2});
                    TN = cat(2,TN,{t_str1,t_str2});
                end
            end
        end
        
        function [H,Y,TN] =bac_fig_53P2()
            method_key = {'factor_IC','factor_MHKQ','S53P2basic','S53P2final'};
            method_info = containers.Map(method_key,...
                {'S53P1IC加权','S53P1MHKQ加权','S53P2全行业IC加权','S53P2分行业IC加权'});
            tN = 'S37.S53_return';
            %var_info = {'method_info','tradingdate','r','r300','r500','dr300','dr500'};
            index_symbol = {'A','000905','000300'};
            index_info = containers.Map(index_symbol,{'A股股票池','500股票池','300股票池'});
            %index result
            sql_str = 'select tradeDate,CHGPct from yuqerdata.yq_index where symbol = "%s" order by tradeDate';
            index_chg = cell(size(index_symbol));
            for i = 1:length(index_symbol)
                if strcmp(index_symbol{i},'A')
                    temp = '000985';
                else
                    temp = index_symbol{i};
                end
                index_chg{i} = fetchmysql(sprintf(sql_str,temp),2);
            end

            sql_str1 = 'select * from %s where mID = "%s" and pID = "%s" order by tradeDate';

            H=[];
            Y={};
            TN = {};
            for m_num = 1:2

                m_id = method_key;
                m_id = m_id{m_num};
                for j = 2
                    p_id=index_symbol{1};
                    p_id1=index_symbol{j};
                    x = fetchmysql(sprintf(sql_str1,tN,m_id,p_id),2);
                    tref = x(:,3);
                    [tref,ia,ib] = intersect(tref,index_chg{j}(:,1));
                    x = cell2mat(x(ia,end));
                    y = cell2mat(index_chg{j}(ib,2));
                    x(1) = 0;
                    y(1) = 0;

                    yc1 = cumprod(1+x);
                    yc2 = cumprod(1+x-y);        
                    t_str = sprintf('%s-%s',method_info(m_id),index_info(p_id1));        
                    t_str1 = sprintf('%s-做多',t_str);
                    t_str2 = sprintf('%s-指数对冲',t_str);        

                    h1 = figure_S53(yc1,tref,t_str1);
                    h2 = figure_S53(yc2,tref,t_str2);

                    H = cat(2,H,[h1,h2]);
                    Y = cat(2,Y,{yc1,yc2});
                    TN = cat(2,TN,{t_str1,t_str2});
                end
            end
        end
        
        
        function [H,re] = bac_figure3_update()
            sql_tmp = 'select * from S37.S51_return where method_info = "S51P3" order by tradingdate';
            y_re = fetchmysql(sql_tmp,2);
            r_s = cell2mat(y_re(:,3));
            tref_s = y_re(:,2);
            tref_s_num = datenum(tref_s);
            p_t = fetchmysql('select distinct(publishDate) from S37.S51_hold_info_p1 where method_info = "S51P3" order by publishDate',2);
            p_t_id = datenum(p_t);

            sql_str = ['select tradeDate,CloseIndex from yuqerdata.yq_index where symbol = "000300" ',...
                'and tradeDate>="%s" order by tradeDate'];
            y_ref = fetchmysql(sprintf(sql_str,'2001-01-01'),2);
            tref_ref = y_ref(:,1);
            tref_ref_num = datenum(tref_ref);
            r_ref = cell2mat(y_ref(:,2));
            r_ref(2:end) = r_ref(2:end)./r_ref(1:end-1)-1;
            r_ref(1) = 0;

            %
            T  = length(p_t);
            wid = 20;
            x = cell(T,1);
            for i = 1:T
                sub_wid = find(tref_ref_num<=p_t_id(i),20,'last');
                sub_r = r_ref(sub_wid);
                sub_r = cumprod(1+sub_r)-1;
                sub_r = sub_r(end);
                if i <T
                    sub_wid2 = tref_s_num>p_t_id(i)&tref_s_num<=p_t_id(i+1);
                else
                    sub_wid2 = tref_s_num>p_t_id(i);
                end
                if sub_r<0.10       
                    sub_re = [tref_s(sub_wid2),num2cell(r_s(sub_wid2))];
                else
                    sub_wid2 = find(sub_wid2);
                    if length(sub_wid2)<=wid
                        sub_re = [tref_s(sub_wid2),num2cell(r_s(sub_wid2))];
                    else
                        sub_wid2_0 = sub_wid2(1:wid);
                        sub_re_0 = [tref_s(sub_wid2_0),num2cell(r_s(sub_wid2_0))];
                        sub_wid2_1 = sub_wid2(wid+1:end);
                        sub_t = tref_s_num(sub_wid2_1);
                        sub_wid2_3 = tref_ref_num>=min(sub_t) & tref_ref_num<=max(sub_t);
                        sub_re_1 = [tref_ref(sub_wid2_3),num2cell(r_ref(sub_wid2_3))];
                        sub_re = [sub_re_0;sub_re_1];
                    end
                end
                x{i} = sub_re';   
            end
            y=[x{:}]';

            %y_ref2=[y_re(:,[2,end]),cumprod(1+cell2mat(y(:,2)))];
            %plot(y_ref2,'LineWidth',2)
            y_r = cell2mat([y_re(:,3:4),y(:,end)]);
            y_r(1,:) = 0;
            y_ref2 = cumprod(1+y_r);
            t_str = cellfun(@(x) strjoin(strsplit(x,'-'),''),y(:,1),'UniformOutput',false);
            T = length(t_str);
            %cut_v = [0,1,3,5];
            leg_str = {'重仓股方法-p3','300指数','重仓股方法择时'};
            H = zeros(2,1);
            H(1) = figure;
            plot(y_ref2,'LineWidth',2);
            set(gca,'xlim',[0,T]);
            set(gca,'XTick',floor(linspace(1,T,15)));
            set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
            set(gca,'XTickLabelRotation',90)    
            setpixelposition(gcf,[223,365,1345,420]);
            legend(leg_str,'Location','best')
            box off

            H(2) = figure;
            temp = [y_r(:,1)-y_r(:,2),y_r(:,3)-y_r(:,2)];
            %temp = temp(:,[1,3]);
            y_re2 = cumprod(1+temp);
            plot(y_re2,'LineWidth',2);
            set(gca,'xlim',[0,T]);
            set(gca,'XTick',floor(linspace(1,T,15)));
            set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
            set(gca,'XTickLabelRotation',90)    
            setpixelposition(gcf,[223,365,1345,420]);
            legend({'p3','p3择时'},'Location','best')
            box off
            y_com = {y_ref2(:,3),y_re2(:,2)};
            t_com = {'重仓股方法择时','重仓股方法择时-对冲300'};
            re = {y_com,t_com};
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

function h = figure_S53(y_re,tref,t_str1)
t_str = cellfun(@(x) strjoin(strsplit(x,'-'),''),tref,'UniformOutput',false);
T = length(t_str);
h = figure;
plot(y_re,'LineWidth',2);
set(gca,'xlim',[0,T]);
set(gca,'XTick',floor(linspace(1,T,15)));
set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
set(gca,'XTickLabelRotation',90)    
setpixelposition(h,[223,365,1345,420]);
%legend(leg_str,'Location','best')
box off
title(t_str1)
end