%{
以下都是万分之几
香港买11卖11
中国买2卖12
美国澳洲日本买1卖1
台湾买2卖32

%}
classdef bac_result_S43<handle
    methods
        function get_all_results(obj)
            obj.update_signal()
            [H1,re1,info1] = obj.bac_re_other();
            [H2,re2,info2] = obj.bac_re_csi();
            key_str = 'S43项目';
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
        function update_signal()
            sprintf('开始执行43升级信号程序')
            dos('python bac_toolS43_V4.py')
        end
        function [H,Y,info] = bac_re_other()
            key_str = sprintf('S43 双底策略指数成分股组合%s',datestr(now,'yyyymmdd'));
            %forex-day,as51,topix
            index_pool = {'hsce','US','HK','forex_day','as51','topix','twse','hk_ggt','kosdaq', ...
                'kospi', 'msci', 'ndx', 'nifty', 'nky', 'RTY', 'set50', 'sx5e',...
                       'ukx', 'xin9i'};
            %index_pool = {'HK'};
            T_index_pool = length(index_pool);
            Y = cell(1,T_index_pool);
            H = zeros(1,T_index_pool);
            info = Y;
            for i = 1:T_index_pool
                index0 =index_pool{i};
                sql_str = sprintf('select ticker,tradeDate,r,sig  from S37.s43_%s order by tradeDate',index0);
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

                t_str = cellfun(@(x) [x(1:4),x(6:7),x(9:10)],tref,'UniformOutput',false);
                T = length(t_str);

                title_str = index0;
                title_str(strfind(title_str,'_')) = '-';

                sub_symbols = symbols;

                [~,ia] = intersect(symbols,sub_symbols);
                sub_Y = X(:,ia);                
                sub_r = mean(sub_Y,2);
                r_c = cumprod(1+sub_r);
                h=figure;
                plot(r_c,'-','LineWidth',2);
                set(gca,'xlim',[0,T]);
                set(gca,'XTick',floor(linspace(1,T,15)));
                set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
                set(gca,'XTickLabelRotation',90)    
                setpixelposition(h,[223,365,1345,420]);
                box off
                title(title_str)    
                sprintf('%s %d-%d',key_str,i,T_index_pool)
                Y{i} = r_c;
                H(i) = h;
                info{i} = title_str;
            end
        end
        function [H,Y,info] = bac_re_csi()
            key_str = sprintf('S43 双底策略指数成分股组合%s',datestr(now,'yyyymmdd'));            
            sql_str = 'select ticker,tradeDate,r,sig  from S37.s43_csi order by tradeDate';
            x0 = fetchmysql(sql_str,2);
            x = x0(:,[1,2,4]);
            x = cell2table(x,'VariableNames',{'s','t','v'});
            pool_w = unstack(x,'v','s');
            pool_w.Properties.RowNames = pool_w.t;
            pool_w.t = [];
            pn0 = write_signal_check();            
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
                r_c = cumprod(1+sub_r);
                h=figure;
                plot(r_c,'-','LineWidth',2);
                set(gca,'xlim',[0,T]);
                set(gca,'XTick',floor(linspace(1,T,15)));
                set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
                set(gca,'XTickLabelRotation',90)    
                setpixelposition(h,[223,365,1345,420]);
                box off
                title(title_str)
                sprintf('%s %d-%d',key_str,i,T_index_pool)
                Y{i} = r_c;
                H(i) = h;
                info{i} = title_str;
            end
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


function pn0 = write_signal_check()
pn0 = fullfile(pwd,'para_pool_S43');
if ~exist(pn0,'dir')
    mkdir(pn0);
end
end