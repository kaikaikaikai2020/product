classdef bac_result_S40<handle
    methods
        function get_all_results(obj)
            obj.M_YCZmin_section_to_series();
            obj.M_Astock_Indicator();
            obj.M_future_CF();
            obj.M_Astock_Future_Indicator();
            obj.M_Foreign_exchange();
            obj.M_Astock_final();
        end
    end
    methods(Static)
        function M_YCZmin_section_to_series()
            keystr = '预测者 截面转时间序列';
            dN= 'ycz_min_series';

            dN_source = 'ycz_min_history';
            symbols = yq_methods.get_symbol_A();
            ycz_symbols = symbols;
            id1 = cellfun(@(x) strcmp(x(1),'6'),symbols);
            ycz_symbols(id1) = cellfun(@(x) ['sh',x],ycz_symbols(id1),'UniformOutput',false);
            ycz_symbols(~id1) = cellfun(@(x) ['sz',x],ycz_symbols(~id1),'UniformOutput',false);

            %获取预测者section数据
            tns_sources = fetchmysql(sprintf('show tables from %s',dN_source),2);
            tns_sources = sort(tns_sources);

            T_tns_source = length(tns_sources);
            sql_str = 'select distinct(symbol) from %s.`%s` order by symbol';
            sql_str_f1 = 'insert into %s.%s select * from %s.`%s`  where symbol = ''%s''';

            if exist('S40_Astock_min_date.mat','file')
                load('S40_Astock_min_date');
                num0 = find(strcmp(tns_sources,sub_tns_sources))+1;
            else
                sub_tns_sources = [];
                num0 = 1;
            end
            for i = num0:T_tns_source
                sub_tns_sources = tns_sources{i};
                sub_symbols = fetchmysql(sprintf(sql_str,dN_source,sub_tns_sources),2);
                T_symbols = length(sub_symbols);
                parfor j = 1:T_symbols
                    %借用mysql语句
                    sub_sql_str = sprintf(sql_str_f1,dN,sub_symbols{j},dN_source,sub_tns_sources,sub_symbols{j});
                    try
                        exemysql(sub_sql_str);
                    catch e_info
                        sprintf(e_info.message)
                    end
                    sprintf('%s : %d-%d  %d -%d',keystr,j,T_symbols,i,T_tns_source)
                end     
            end
            if ~isempty(sub_tns_sources)
                save('S40_Astock_min_date.mat','sub_tns_sources');
            end
            
        end
        
        function M_Astock_Indicator()
            key_str = 'S40SMT策略指数回测';
            symbols = {'sh000016','sh000905','sh000300'};
            symbols_index = {'50','500','300'};
            T_symbols = length(symbols);
            sta_re = cell(T_symbols,1);
            error_ind = zeros(T_symbols,1);
            write_sel = true;
            if write_sel
                pn_write = fullfile(pwd,'计算结果');
                if ~exist(pn_write,'dir')
                    mkdir(pn_write)
                end
                obj_wd = wordcom(fullfile(pn_write,sprintf('%s.doc',key_str)));
                xls_fn = fullfile(pn_write,sprintf('%s.xlsx',key_str));
            end

            for i_sym = 1:T_symbols
                %参数设置
                try 
                    P = [];
                    P.feeOpen=5/100000;
                    P.feeClose=5/100000;
                    P.matchRecord=1;%匹配数据源：沪深300
                    P.tradeRecord=1;%交易数据源：股指期货主力合约
                    P.tradeMin=120;%使用早盘120分钟K线数据进行分形匹配
                    P.dayMin=240;%每个交易日共240根1分钟K线
                    P.M=20;%找M个最为相似的交易日
                    P.muchPara=0.5;%多数上涨或下跌比例
                    P.deanMethod=3;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离
                    P.stopMethod=3;%1收盘价止损/2触发价止损，不考虑能否交易/3跳开则开盘价止损，否则触发价止损
                    P.testStart='2010-4-16';
                    P.trade_mode = 2;%1只多仓 2 多仓和空仓

                    title_str = symbols{i_sym};
                    sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
                        'hour(tradingdate),minute(tradingdate),open,close from S40.%s'];
                    sub_sql_str = sprintf(sql_str,title_str);
                    x = fetchmysql(sub_sql_str);
                    if ~isempty(x)
                        sub_tt = datenum([x(end,1:5),0]);
                        sub_tt = datestr(sub_tt,'yyyy-mm-dd HH:MM:SS');
                        sql_str1 = ['select t_year,t_month,t_day,t_hour,t_minute,open,close from ',...
                            'pytdx_data.tdx_min_%s where tradingdate >''%s'' order by tradingdate'];%tdx_min_
                        sub_sql_str1 = sprintf(sql_str1,symbols_index{i_sym},sub_tt);
                        x1 = fetchmysql(sub_sql_str1);
                        x = cat(1,x,x1);
                    else
                        sql_str1 = ['select t_year,t_month,t_day,t_hour,t_minute,open,close from ',...
                            'pytdx_data.tdx_min_%s order by tradingdate'];%tdx_min_
                        sub_sql_str1 = sprintf(sql_str1,symbols_index{i_sym});
                        x = fetchmysql(sub_sql_str1);
                    end

                    temp = x(:,4)*100+x(:,5);
                    id_sel = temp>930&temp<=1500;
                    x = x(id_sel,:);
                    if isempty(x)
                        sprintf('%s %s 时间不足6年，跳出',key_str,title_str)
                        continue
                    end
                    t = datenum([x(:,1:5),zeros(size(x(:,1)))]);
                    %统计中间有停牌的情况，并剔除
                    day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
                    day_tick_u = unique(day_tick);
                    ind_miss = false(size(day_tick));
                    ind_miss_u = false(size(day_tick_u));
                    T = length(day_tick_u);
                    for i = 1:T
                        sub_ind = eq(day_tick,day_tick_u(i));
                        if sum(sub_ind)<240
                            ind_miss(sub_ind) = true;
                            ind_miss_u(i) = true;
                        end
                    end

                    day_tick_u(ind_miss_u) = [];
                    x(ind_miss,:) = [];
                    t(ind_miss,:) = [];

                    %初始时间设定
                    min_day_num = 210*6;
                    if length(day_tick_u)<min_day_num
                        sprintf('%s %s 时间不足6年，跳出',key_str,title_str)
                        continue
                    else
                        temp = num2str(day_tick_u(min_day_num/2+1));
                        P.testStart=sprintf('%s-%s-%s',temp(1:4),temp(5:6),temp(7:8));
                    end

                    openprice = x(:,end-1);
                    closeprice = x(:,end);
                    %SMTTradingModelTool(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate)
                    [tradeYield,result1,tradeDetail,yearDetail,h] = SMTTradingModelTool2(closeprice,openprice,t,closeprice,openprice,t,P);
                    ah=gca;
                    title(ah,title_str)
                    if write_sel
                        obj_wd.pasteFigure(h,title_str);
                    end
                    y_c = cumprod(tradeYield(:,2)+1);
                    %统计参数
                    [v0,v_str0] = curve_static(y_c,[],false);
                    [v,v_str] = ad_trans_sta_info(v0,v_str0); 
                    result2 = [v_str;v]';
                    result = [{'',title_str};[result1;result2]];
                    sta_re{i_sym} = result;
                    sprintf('%s %d-%d',key_str,i_sym,T_symbols)
                catch
                    error_ind(i_sym) = 1;
                    sprintf('Error %s %d-%d',key_str,i_sym,T_symbols)
                end

            end
            y = [sta_re{:}];
            y = y(:,[1,2:2:end]);
            y = y';
            if write_sel
                obj_wd.CloseWord();
                xlswrite(xls_fn,y);
            end
        end
        
        function M_future_CF()
            info = ['select secFullName,contractObject,exchangeCD  ',...
                'from yuqerdata.yq_FutuGet group by contractObject order by contractObject'];
            info = fetchmysql(info,2);

            key_str = 'S40SMT策略国内商品期回测结果';

            dN = 'Future_min_Jinshuyuan';
            dN_tdx = 'future_min_data';
            symbols = fetchmysql(sprintf('show tables from %s',dN),2);
            T_symbols = length(symbols);
            sta_re = cell(T_symbols,1);
            error_ind = zeros(T_symbols,1);

            write_sel = true;
            if write_sel
                pn_write = fullfile(pwd,'计算结果');
                if ~exist(pn_write,'dir')
                    mkdir(pn_write)
                end
                obj_wd = wordcom(fullfile(pn_write,sprintf('%s.doc',key_str)));
                xls_fn = fullfile(pn_write,sprintf('%s.xlsx',key_str));
            end

            tradeDay_av = zeros(T_symbols,1);
            for i_sym = 1:T_symbols

                P = [];
                P.feeOpen=5/100000;
                P.feeClose=5/100000;
                P.matchRecord=1;%匹配数据源：沪深300
                P.tradeRecord=1;%交易数据源：股指期货主力合约
                P.tradeMin=135;%使用早盘135分钟K线数据进行分形匹配
                P.dayMin=225;%每个交易日共225根1分钟K线
                P.M=20;%找M个最为相似的交易日
                P.muchPara=0.5;%多数上涨或下跌比例
                P.deanMethod=3;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离
                P.stopMethod=3;%1收盘价止损/2触发价止损，不考虑能否交易/3跳开则开盘价止损，否则触发价止损
                P.testStart='2010-4-16';
                P.trade_mode = 2;%1只多仓 2 多仓和空仓

                title_str = symbols{i_sym};
                temp_id = strcmpi(info(:,2),title_str);
                if any(temp_id)
                    title_str2 = info{temp_id,1};
                    tn_tdx = sprintf('%s_%s',info{temp_id,3},info{temp_id,2});
                else
                    title_str2 = title_str;
                    tn_tdx = [];
                end

                sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
                    'hour(tradingdate),minute(tradingdate),openprice,closeprice from %s.%s ',...
                    'order by tradingdate'];
                sub_sql_str = sprintf(sql_str,dN,title_str);
                x = fetchmysql(sub_sql_str);
                if ~isempty(tn_tdx)
                    sub_t_max = x(end,1:6);
                    sub_t_max(end) = 0;
                    sub_t_max = datestr(datenum(sub_t_max),'yyyy-mm-dd HH:MM:SS');
                    sql_str_tdx = ['select t_year,t_month,t_day,',...
                        't_hour,t_minute,open,close from %s.%s ',...
                        'where tradingdate>''%s'' order by tradingdate'];
                    sub_sql_str = sprintf(sql_str_tdx,dN_tdx,tn_tdx,sub_t_max);
                    x1 = fetchmysql(sub_sql_str);
                    x = cat(1,x,x1);
                else
                    sql_str_tdx = ['select t_year,t_month,t_day,',...
                        't_hour,t_minute,open,close from %s.%s ',...
                        ' order by tradingdate'];
                    sub_sql_str = sprintf(sql_str_tdx,dN_tdx,tn_tdx);
                    x = fetchmysql(sub_sql_str);
                end
                if isempty(x)
                    sprintf('%s %s 时间不足6年，跳出',key_str,title_str)
                    continue
                end
                temp = x(:,4)*100+x(:,5);
                id_sel = temp>900&temp<=1500;
                x = x(id_sel,:);
                if isempty(x)
                    continue
                end
                t = datenum([x(:,1:5),zeros(size(x(:,1)))]);
                %统计中间有停牌的情况，并剔除
                day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
                day_tick_u = unique(day_tick);
                ind_miss = false(size(day_tick));
                ind_miss_u = false(size(day_tick_u));
                T = length(day_tick_u);
                temp_tradedays = zeros(T,1);
                for i = 1:T
                    sub_ind = eq(day_tick,day_tick_u(i));
                    if sum(sub_ind)<225
                        ind_miss(sub_ind) = true;
                        ind_miss_u(i) = true;
                    end
                    temp_tradedays(i) = sum(sub_ind);
                end
                temp_tradedays(ind_miss_u) = [];
                if ~isempty(temp_tradedays)
                    tradeDay_av(i_sym) = mean(temp_tradedays);
                end

                day_tick_u(ind_miss_u) = [];
                x(ind_miss,:) = [];
                t(ind_miss,:) = [];

                %初始时间设定
                min_day_num = 210*6;
                if length(day_tick_u)<min_day_num
                    sprintf('%s %s 时间不足6年，跳出',key_str,title_str)
                    continue
                else
                    temp = num2str(day_tick_u(min_day_num/2+1));
                    P.testStart=sprintf('%s-%s-%s',temp(1:4),temp(5:6),temp(7:8));
                end

                openprice = x(:,end-1);
                closeprice = x(:,end);
                %SMTTradingModelTool(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate)
                [tradeYield,result1,tradeDetail,yearDetail,h] = SMTTradingModelTool2(closeprice,openprice,t,closeprice,openprice,t,P);
                ah = gca;
                title(ah,title_str2);
                if write_sel
                    obj_wd.pasteFigure(h,title_str2);
                end
                y_c = cumprod(tradeYield(:,2)+1);
                %统计参数
                [v0,v_str0] = curve_static(y_c,[],false);
                [v,v_str] = ad_trans_sta_info(v0,v_str0); 
                result2 = [v_str;v]';
                result = [{'',title_str2};[result1;result2]];
                sta_re{i_sym} = result;
                sprintf('%s %d-%d',key_str,i_sym,T_symbols)

            end
            y = [sta_re{:}];
            y = y(:,[1,2:2:end]);
            y = y';
            if write_sel
                obj_wd.CloseWord()
                xlswrite(xls_fn,y);
            end
        end
        
        function M_Astock_Future_Indicator()
            key_str = 'S40SMT策略股指期货测试';
            symbols = {'ZJIH','ZJIC','ZJIF'};
            symbols_tdx = {'ccfx_ih','ccfx_ic','ccfx_if'};
            dN_tdx = 'Future_min_data';
            T_symbols = length(symbols);
            sta_re = cell(T_symbols,1);
            error_ind = zeros(T_symbols,1);

            write_sel = true;
            if write_sel
                pn_write = fullfile(pwd,'计算结果');
                if ~exist(pn_write,'dir')
                    mkdir(pn_write)
                end
                obj_wd = wordcom(fullfile(pn_write,sprintf('%s.doc',key_str)));
                xls_fn = fullfile(pn_write,sprintf('%s.xlsx',key_str));
            end
            for i_sym = 1:T_symbols
                %参数设置
                P = [];
                P.feeOpen=5/100000;
                P.feeClose=5/100000;
                P.matchRecord=1;%匹配数据源：沪深300
                P.tradeRecord=1;%交易数据源：股指期货主力合约
                P.tradeMin=120;%使用早盘120分钟K线数据进行分形匹配
                P.dayMin=240;%每个交易日共240根1分钟K线
                P.M=20;%找M个最为相似的交易日
                P.muchPara=0.5;%多数上涨或下跌比例
                P.deanMethod=3;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离
                P.stopMethod=3;%1收盘价止损/2触发价止损，不考虑能否交易/3跳开则开盘价止损，否则触发价止损
                P.testStart='2010-4-16';
                P.trade_mode = 2;%1只多仓 2 多仓和空仓

                title_str = symbols{i_sym};
                sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
                    'hour(tradingdate),minute(tradingdate),open,close from S40.%s'];
                sub_sql_str = sprintf(sql_str,title_str);
                x = fetchmysql(sub_sql_str);
                if isempty(x)
                    sql_str_tdx = ['select t_year,t_month,t_day,',...
                        't_hour,t_minute,open,close from %s.%s ',...
                        ' order by tradingdate'];
                    sub_sql_str = sprintf(sql_str_tdx,dN_tdx,symbols_tdx{i_sym});
                    x = fetchmysql(sub_sql_str);
                else
                    sub_t_max = x(end,1:6);
                    sub_t_max(end) = 0;
                    sub_t_max = datestr(datenum(sub_t_max),'yyyy-mm-dd HH:MM:SS');
                    sql_str_tdx = ['select t_year,t_month,t_day,',...
                        't_hour,t_minute,open,close from %s.%s ',...
                        'where tradingdate>''%s'' order by tradingdate'];
                    sub_sql_str = sprintf(sql_str_tdx,dN_tdx,symbols_tdx{i_sym},sub_t_max);
                    x1 = fetchmysql(sub_sql_str);
                    x = cat(1,x,x1);
                end
                if isempty(x)
                    sprintf('%s %s 时间不足6年，跳出',key_str,title_str)
                    continue
                end
                temp = x(:,4)*100+x(:,5);
                id_sel = temp>930&temp<=1500;
                x = x(id_sel,:);

                t = datenum([x(:,1:5),zeros(size(x(:,1)))]);
                %统计中间有停牌的情况，并剔除
                day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
                day_tick_u = unique(day_tick);
                ind_miss = false(size(day_tick));
                ind_miss_u = false(size(day_tick_u));
                T = length(day_tick_u);
                for i = 1:T
                    sub_ind = eq(day_tick,day_tick_u(i));
                    if sum(sub_ind)<240
                        ind_miss(sub_ind) = true;
                        ind_miss_u(i) = true;
                    end
                end

                day_tick_u(ind_miss_u) = [];
                x(ind_miss,:) = [];
                t(ind_miss,:) = [];

                %初始时间设定
                min_day_num = 210*4;
                if length(day_tick_u)<min_day_num
                    sprintf('%s %s 时间不足6年，跳出',key_str,title_str)
                    continue
                else
                    temp = num2str(day_tick_u(min_day_num/2+1));
                    P.testStart=sprintf('%s-%s-%s',temp(1:4),temp(5:6),temp(7:8));
                end

                openprice = x(:,end-1);
                closeprice = x(:,end);
                %SMTTradingModelTool(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate)
                [tradeYield,result1,tradeDetail,yearDetail,h] = SMTTradingModelTool2(closeprice,openprice,t,closeprice,openprice,t,P);
                ah=gca;
                title(ah,title_str)
                if write_sel
                    obj_wd.pasteFigure(h,title_str);
                end
                y_c = cumprod(tradeYield(:,2)+1);
                %统计参数
                [v0,v_str0] = curve_static(y_c,[],false);
                [v,v_str] = ad_trans_sta_info(v0,v_str0); 
                result2 = [v_str;v]';
                result = [{'',title_str};[result1;result2]];
                sta_re{i_sym} = result;
                sprintf('%s %d-%d',key_str,i_sym,T_symbols)   
            end
            y = [sta_re{:}];
            y = y(:,[1,2:2:end]);
            y = y';
            if write_sel
                obj_wd.CloseWord()
                xlswrite(xls_fn,y);
            end
        end
        
        function M_Foreign_exchange()
            write_sel = true;
            key_str = 'S40SMT策略主要外汇';

            dN = 'foreign_index_min_V2';
            symbols = strsplit('cadchf,cadjpy,chfjpy,eurcad,eurchf,eurgbp,eurjpy,eurusd,gbpcad,gbpchf,gbpjpy,usdcad,usdchf,audchf,audjpy,euraud,eurnzd',',');
            symbols1 = strsplit('audcad,audchf,audjpy,audnzd,audusd,nzdusd,nzdjpy,euraud,eurnzd,usdjpy',',');
            T_symbols = length(symbols);
            re = cell(T_symbols,1);
            X = re;
            parfor i_sym =  1:T_symbols

                P = [];
                P.feeOpen=1/100000;
                P.feeClose=1/100000;
                P.matchRecord=1;%匹配数据源：沪深300
                P.tradeRecord=1;%交易数据源：股指期货主力合约
                P.tradeMin=360;%使用早盘135分钟K线数据进行分形匹配
                P.dayMin=360;%每个交易日共225根1分钟K线
                P.M=20;%找M个最为相似的交易日
                P.muchPara=0.5;%多数上涨或下跌比例
                P.deanMethod=3;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离
                P.stopMethod=3;%1收盘价止损/2触发价止损，不考虑能否交易/3跳开则开盘价止损，否则触发价止损
                P.testStart='2010-4-16';
                P.trade_mode = 2;%1只多仓 2 多仓和空仓
                P.out_sel = false;

                title_str = symbols{i_sym};  
                if any(strcmp(symbols1,title_str))
                    t_num1 = 200;
                    t_num2 = 1400;
                else
                    t_num1 = 800;
                    t_num2 = 2000;
                end
                sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
                    'hour(tradingdate),minute(tradingdate),open,close from %s.%s ',...
                    ' order by tradingdate'];
                sub_sql_str = sprintf(sql_str,dN,title_str);
                x = fetchmysql(sub_sql_str);
                if isempty(x)
                    sprintf('%s %s 时间不足4年，跳出',key_str,title_str)
                    continue
                end
                temp = x(:,4)*100+x(:,5);
                id_sel = temp>t_num1&temp<=t_num2;
                x = x(id_sel,:);
                if isempty(x)
                    continue
                end
                t = datenum([x(:,1:5),zeros(size(x(:,1)))]);
                %统计中间有停牌的情况，并剔除
                day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
                day_tick_u = unique(day_tick);
                ind_miss = false(size(day_tick));
                ind_miss_u = false(size(day_tick_u));
                T = length(day_tick_u);
                temp_tradedays = zeros(T,1);
                for i = 1:T
                    sub_ind = eq(day_tick,day_tick_u(i));
                    if ~eq(sum(sub_ind),720)
                        ind_miss(sub_ind) = true;
                        ind_miss_u(i) = true;
                    end
                    temp_tradedays(i) = sum(sub_ind);
                end
                temp_tradedays(ind_miss_u) = [];

                day_tick_u(ind_miss_u) = [];
                x(ind_miss,:) = [];
                t(ind_miss,:) = [];

                %初始时间设定
                min_day_num = 210*4;
                if length(day_tick_u)<min_day_num
                    sprintf('%s %s 时间不足4年，跳出',key_str,title_str)
                    continue
                else
                    temp = num2str(day_tick_u(min_day_num/2+1));
                    P.testStart=sprintf('%s-%s-%s',temp(1:4),temp(5:6),temp(7:8));
                end

                openprice = x(:,end-1);
                closeprice = x(:,end);
                [tradeYield,result1,tradeDetail,yearDetail,h,assure_ratio] = SMTTradingModelTool2(closeprice,openprice,t,closeprice,openprice,t,P);
                re{i_sym} = [tradeYield,assure_ratio];
                X{i_sym} = tradeYield(:,1:2)';
                sprintf('%s %d-%d',key_str,i_sym,T_symbols)

            end

            t = [X{:}]';
            t = unique(t(:,1));

            r = zeros(length(t),T_symbols);
            for i = 1:T_symbols
                temp = X{i}';
                [~,ia,ib] = intersect(t,temp(:,1));
                r(ia,i) = temp(ib,2);
            end
            r_c = cumprod(1+mean(r,2));
            title_str = '组合';
            t_str = cellstr(datestr(t,'yyyymmdd'));
            T = length(t_str);
            h=figure;
            subplot(2,1,1)
            plot(r_c,'-','LineWidth',2);
            set(gca,'xlim',[0,T]);
            set(gca,'XTick',floor(linspace(1,T,15)));
            set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
            set(gca,'XTickLabelRotation',90)    
            setpixelposition(h,[223,365,1345,600]);
            box off
            title(sprintf('%s-组合曲线',title_str))
            %每年收益
            t_y = year(t);
            t_y_u = unique(t_y);
            r_year = zeros(size(t_y_u));
            for j = 1:length(t_y_u)
            sub_r = r(eq(t_y,t_y_u(j)),:);
            temp = cumprod(1+mean(sub_r,2));
            r_year(j) = temp(end)-1;
            end
            subplot(2,1,2)
            bar(t_y_u,r_year)
            box off
            title(sprintf('%s-每年收益统计',title_str))

            %三条曲线的参数
            r_1 = [r,mean(r,2)];
            r_str = [symbols,'组合'];
            sub_re = cell(T_symbols,1);
            for j = 1:T_symbols+1
                sub_y = cumprod(1+r_1(:,j));
                [v0,v_str0] = curve_static(sub_y,[],false);
                [v,v_str] = ad_trans_sta_info(v0,v_str0);
                if eq(j,1)
                    sub_re{j} = [[{''},r_str(j)];[v_str;v]'];
                else
                    sub_re{j} = [r_str(j);v'];
                end
            end
            sub_re = [sub_re{:}]';

            if write_sel
                pn_write = fullfile(pwd,'计算结果');
                if ~exist(pn_write,'dir')
                    mkdir(pn_write)
                end
                obj_wd = wordcom(fullfile(pn_write,sprintf('%s.doc',key_str)));
                xls_fn = fullfile(pn_write,sprintf('%s.xlsx',key_str));

                obj_wd.pasteFigure(h,key_str);
                obj_wd.CloseWord();
                xlswrite(xls_fn,sub_re);
            end
            
        end
        
        function M_Astock_final()
            tn = 'S40.A_stock_signal';
            var_info={'symbol', 'tradenum', 'f_val'};

            key_str = 'S40SMT策略A股计算';
            symbols = fetchmysql('show tables from ycz_min_series',2);
            T_symbols = length(symbols);
            sta_re = cell(T_symbols,1);
            sta_re2 = sta_re;
            Y1 = sta_re;
            Y2 = sta_re;
            error_ind = zeros(T_symbols,1);
            %sql_str_r = 'select tradeDate,chgPct from yuqerdata.yq_dayprice where symbol = ''%s'' order by tradeDate';
            sql_str_r = 'select tradeDate,closeprice/openprice-1 from yuqerdata.yq_dayprice where symbol = ''%s'' order by tradeDate';
            X = cell(T_symbols,1);
            parfor i_sym =1:T_symbols
                %参数设置
                P = [];
                P.feeOpen=1.5/1000/2;
                P.feeClose=1.5/1000/2;
                P.matchRecord=1;%匹配数据源：沪深300
                P.tradeRecord=1;%交易数据源：股指期货主力合约
                P.tradeMin=240;%使用早盘120分钟K线数据进行分形匹配
                P.dayMin=240;%每个交易日共240根1分钟K线
                P.M=20;%找M个最为相似的交易日
                P.muchPara=0.5;%多数上涨或下跌比例
                P.deanMethod=3;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离/5 dynamic time warping
                P.stopMethod=3;%1收盘价止损/2触发价止损，不考虑能否交易/3跳开则开盘价止损，否则触发价止损
                P.testStart='2010-4-16';
                P.trade_mode = 1;%1只多仓 2 多仓和空仓
                P.cut_return = -inf;

                title_str = symbols{i_sym};
                sql_str = ['select year(tradingdate),month(tradingdate),day(tradingdate),',...
                    'hour(tradingdate),minute(tradingdate),open,close from ycz_min_series.%s'];
                sub_sql_str = sprintf(sql_str,title_str);
                x = fetchmysql(sub_sql_str);
                if isempty(x)
                    sprintf('%s %s 时间不足6年，跳出',key_str,title_str)
                    continue
                end
                t = datenum([x(:,1:5),zeros(size(x(:,1)))]);
                %统计中间有停牌的情况，并剔除
                day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
                day_tick_u = unique(day_tick);
                ind_miss = false(size(day_tick));
                ind_miss_u = false(size(day_tick_u));
                T = length(day_tick_u);
                for i = 1:T
                    sub_ind = eq(day_tick,day_tick_u(i));
                    if sum(sub_ind)<240
                        ind_miss(sub_ind) = true;
                        ind_miss_u(i) = true;
                    end
                end

                day_tick_u(ind_miss_u) = [];
                x(ind_miss,:) = [];
                t(ind_miss,:) = [];

                %初始时间设定
                min_day_num = 210*6;
                if length(day_tick_u)<min_day_num
                    sprintf('%s %s 时间不足6年，跳出',key_str,title_str)
                    continue
                else
                    temp = num2str(day_tick_u(min_day_num/2+1));

                    sub_sql_str =[ 'select tradenum,f_val from S40.A_stock_signal where ',...
                        'symbol = ''%s'' order by tradenum'];
                    ind_before = fetchmysql(sprintf(sub_sql_str,symbols{i_sym}));
                    if isempty(ind_before)        
                        P.testStart=sprintf('%s-%s-%s',temp(1:4),temp(5:6),temp(7:8));
                    else
                        P.testStart=datestr(ind_before(end,1)+1,'yyyy-mm-dd');
                    end
                end

                openprice = x(:,end-1);
                closeprice = x(:,end);
                %SMTTradingModelTool(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate)
                if datenum(P.testStart)<=t(end)
                    ind0 = SMTTradingModelTool_update1(closeprice,openprice,t,closeprice,openprice,t,P);
                else
                    ind0 = 0;
                end
                ind = [ind_before;ind0];
                %mysql_data
                ind0 = num2cell(ind0(:,[1,1:end]));
                ind0(:,1) = symbols(i_sym);
                X{i_sym} = ind0';
                r = fetchmysql(sprintf(sql_str_r,title_str(3:end)),2);                
                t1= cellstr(datestr(ind(:,1),'yyyy-mm-dd'));
                [sub_tref,ia,ib] = intersect(r(:,1),t1(:,1));
                x=[ind(ib,:),cell2mat(r(ia,end))];

                ind1 = x(:,2);
                ind1(2:end) = ind1(1:end-1);
                %ind1 = -ind1;
                ind1(ind1<0) = 0;
                r = x(:,end);
                r2 = r;
                fee_ind = find(~eq(diff(ind1),0))+1;
                r2(fee_ind) = r2(fee_ind) -P.feeOpen-P.feeClose;
                y_c2 = cumprod(1+r.*ind1);
                y_c = cumprod(1+r2.*ind1);
                %统计参数
                [v0,v_str0] = curve_static(y_c,[],false);
                [v,v_str] = ad_trans_sta_info(v0,v_str0); 
                result2 = [v_str;v]';
                result = [{'',title_str};result2];
                sta_re{i_sym} = result;
                sprintf('%s %d-%d',key_str,i_sym,T_symbols)
                Y1{i_sym} = [sub_tref,num2cell([y_c,y_c2])];
                %反向结果
                ind1 = x(:,2);
                ind1(2:end) = ind1(1:end-1);
                ind1 = -ind1;
                ind1(ind1<0) = 0;
                r = x(:,end);
                r2 = r;
                fee_ind = find(~eq(diff(ind1),0))+1;
                r2(fee_ind) = r2(fee_ind) -P.feeOpen-P.feeClose;
                y_c2 = cumprod(1+r.*ind1);
                y_c = cumprod(1+r2.*ind1);
                %统计参数
                [v0,v_str0] = curve_static(y_c,[],false);
                [v,v_str] = ad_trans_sta_info(v0,v_str0); 
                result2 = [v_str;v]';
                result = [{'',title_str};result2];
                sta_re2{i_sym} = result;
                sprintf('%s %d-%d',key_str,i_sym,T_symbols)
                Y2{i_sym} = [sub_tref,num2cell([y_c,y_c2])];          
            end

            X1 = [X{:}]';
            datainsert_adair(tn,var_info,X1);

            y = [sta_re{:}];
            y = y(:,[1,2:2:end]);
            y = y';

            y2 = [sta_re2{:}];
            y2 = y2(:,[1,2:2:end]);
            y2 = y2';

            write_sel = true;
            if write_sel
                pn_write = fullfile(pwd,'计算结果');
                if ~exist(pn_write,'dir')
                    mkdir(pn_write)
                end
                obj_wd = wordcom(fullfile(pn_write,sprintf('%s.doc',key_str)));
                xls_fn = fullfile(pn_write,sprintf('%s.xlsx',key_str));
            end

            T = length(Y1);
            tref = yq_methods.get_tradingdate();
            Y = zeros(length(tref),T);
            Y_0 = Y;
            for i = 1:T
                sub_x = Y1{i};
                if ~isempty(sub_x)
                    [~,ia,ib] = intersect(tref,sub_x(:,1));
                    temp = cell2mat(sub_x(:,2:3));
                    temp(2:end,:) = temp(2:end,:)./temp(1:end-1,:)-1;
                    temp(1,:) = 0;
                    Y(ia,i) = temp(ib,1);
                    Y_0(ia,i) = temp(ib,2);
                end        
                sprintf('%d-%d',i,T)
            end

            sel_ind1 = sum(abs(Y),1)>1;
            symbols2 = symbols(sel_ind1);
            Y = Y(:,sel_ind1);
            Y_0 = Y_0(:,sel_ind1);
            sel_ind2 = sum(abs(Y),2)>1;
            tref2 = tref(sel_ind2);
            Y = Y(sel_ind2,:);
            Y_0 = Y_0(sel_ind2,:);

            symbols2 = cellfun(@(x) x(3:end),symbols2,'UniformOutput',false);

            t2 = datenum(tref2);
            index_pool = {'000300','000905','000016','000001'};
            index_name = {'沪深300','中证500','上证50','上证综指'};
            sql_str1 = 'select tradeDate,closeIndex from   yuqerdata.yq_index where symbol = ''%s'' order by tradeDate';
            sta_re = cell(size(index_name));
            for index_sel = 1:length(index_pool)
                title_str = index_name{index_sel};
                sub_x = fetchmysql(sprintf(sql_str1,index_pool{index_sel}),2);
                sub_symbols = yq_methods.get_index_pool(index_pool{index_sel},sub_x{end,1});

                %限制股票池
                [~,ia] = intersect(symbols2,sub_symbols);
                [sub_tref,ia1,ib1] = intersect(tref2,sub_x(:,1));

                sub_Y = Y(ia1,ia);
                sub_Y0 = Y_0(ia1,ia);
                sub_t_num = t2(ia1);
                sub_x = sub_x(ib1,:);
                sub_x_c = cell2mat(sub_x(:,2));
                sub_x_c = sub_x_c./sub_x_c(1);

                t_str = sub_tref;
                T = length(t_str);

                r  = mean(sub_Y,2);
                r0 = mean(sub_Y0,2);
                r(1) = 0;
                r0(1) = 0;
                r_c = cumprod(1+r);
                r_c0 = cumprod(1+r0);

                h1=figure;
                subplot(2,1,1)
                plot([r_c,r_c0,sub_x_c],'-','LineWidth',2);
                set(gca,'xlim',[0,T]);
                set(gca,'XTick',floor(linspace(1,T,15)));
                set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
                set(gca,'XTickLabelRotation',90)    
                setpixelposition(h1,[223,365,1345,600]);
                legend({'组合曲线有手续费','组合曲线无手续费',index_pool{index_sel}},'NumColumns',3,'Location','best')
                box off
                title(sprintf('%s-组合曲线',title_str))
                %每年收益
                t_y = year(sub_t_num);
                t_y_u = unique(t_y);
                r_year = zeros(size(t_y_u));
                for j = 1:length(t_y_u)
                sub_r = r(eq(t_y,t_y_u(j)),:);
                temp = cumprod(1+mean(sub_r,2));
                r_year(j) = temp(end)-1;
                end
                subplot(2,1,2)
                bar(t_y_u,r_year)
                box off
                title(sprintf('%s-每年收益统计',title_str))
                if write_sel
                    obj_wd.pasteFigure(h1,title_str);
                end
                %统计参数
                sub_re = cell(2,1);
                sub_y = [sub_x_c,r_c];
                sub_title_str = {title_str,sprintf('%s组合',title_str)};
                for j = 1:2
                    [v0,v_str0] = curve_static(sub_y(:,j),[],false);
                    [v,v_str] = ad_trans_sta_info(v0,v_str0);
                    if eq(j,1)
                        sub_re{j} = [[{''},sub_title_str(j)];[v_str;v]'];
                    else
                        sub_re{j} = [sub_title_str(j);v'];
                    end
                end
                sta_re{index_sel} = [sub_re{:}];
            end
            sta_re = [sta_re{:}]';
            if write_sel
                obj_wd.CloseWord()
                xlswrite(xls_fn,sta_re);
            end
        end
    end
    
    
end