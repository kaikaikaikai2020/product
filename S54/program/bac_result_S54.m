classdef bac_result_S54 < handle
    methods
        function obj = bac_result_S54()
            addpath(genpath(fullfile(pwd,'jplv7')))
            obj.iniA41();
        end
        function [H,RE] = get_all_results(obj)
            tic
            H = {};
            RE = {};
            [h,re] = obj.A33();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A38_EWAEWC();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A38_forex();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A38cis();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A41_bog();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A41_csi();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A43_indexArb();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A44_andrewlo();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A44_csi();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A51_audcad();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A63();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A71();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            [h,re] = obj.A57();
            [H,RE] = app_data_S54(obj,h,re,H,RE);
            T = length(RE);
            for i = 1:T
                sub_X = RE{i};
                if size(sub_X,2)>1
                    RE{i} = sub_X';
                end
            end
            toc
            key_str = 'S54 策略验证';
            file_name = sprintf('%s%s',key_str,datestr(now,'yyyy-mm-dd'));
            pn0 = fullfile(pwd,'计算结果');
            if ~exist(pn0,'dir')
                mkdir(pn0);
            end
            obj_wd = wordcom(fullfile(pn0,sprintf('%s.doc',file_name)));
            H=[H{:}];
            for i = 1:length(H)
                obj_wd.pasteFigure(H(i),' ');
            end
            obj_wd.CloseWord();
            %yc = [re1,re2];
            %yt = [info1,info2];
            temp = [RE{:}]';
            yc = temp(:,1);
            yt = temp(:,2);
            for i = 1:length(yc)
                if eq(yc{i}(1),0)
                    yc{i} = yc{i}+1;
                end
            end
            %yt = [yt{:}];
            sta_re = curve_static_batch(yc,yt);
            xlstocsv_adair(fullfile(pn0,sprintf('%s.xlsx',file_name)),sta_re) 
        end
        function [H,RE] = app_data_S54(obj,h,re,H,RE)
            h = obj.trans_h_size(h);
            H = cat(1,H,{h});
            RE = cat(1,RE,{re'});
        end
    end
    methods(Static)
        function h = trans_h_size(h)
            if size(h,1)>size(h,2)
                h = h';
            end
        end
        function iniA41()            
            sprintf('开始升级A41数据')
            load('inputDataOHLCDaily_stocks_20120424', 'stocks');
            a = sprintf('"%s"',strjoin(stocks,'","'));
            sql_str = ['select ticker,tradeDate,openPrice,highPrice,lowPrice,closePrice from polygon.usastock_day ',...
            'where tradeDate> (select max(tradeDate) from S54.A41data) and ticker in (%s)'];
            x = fetchmysql(sprintf(sql_str,a),2);
            %dN = 'S54';
            tN = 'S54.A41data';
            var_info = {'ticker','tradeDate','openPrice','highPrice','lowPrice','closePrice'};
            if ~isempty(x)
                datainsert_adair(tN,var_info,x)
            end
            sprintf('完成升级A41数据')            
        end
        %A33
        function [h,bac_re] = A33()
            key_str = 'S54-A33';
            sprintf('开始%s',key_str)
            sql_str = 'select tradeDate,closePrice from polygon.usastock_day where ticker = "%s" order by tradeDate';
            x=fetchmysql(sprintf(sql_str,'GLD'),2);
            y=fetchmysql(sprintf(sql_str,'USO'),2);
            [tref,ia,ib] = intersect(x(:,1),y(:,1));
            x = cell2mat(x(ia,2));
            y = cell2mat(y(ib,2));
            lookback=20; % Lookback set arbitrarily short
            hedgeRatio=NaN(size(x, 1), 1);
            for t=lookback:size(hedgeRatio, 1)
                %regression_result=ols(y(t-lookback+1:t), [x(t-lookback+1:t) ones(lookback, 1)]);    
                %hedgeRatio(t)=regression_result.beta(1);
                regression_result=regress(y(t-lookback+1:t), [x(t-lookback+1:t) ones(lookback, 1)]);   
                hedgeRatio(t)=regression_result(1);
            end
            y2=[x y];
            yport=sum([-hedgeRatio ones(size(hedgeRatio))].*y2, 2); % The net market value of the portfolio is same as the "spread"
            hedgeRatio(1:lookback)=[]; % Removed because hedge ratio is indterminate

            if ~isempty(tref)
                tref(1:lookback) = [];
            end

            yport(1:lookback)=[]; 
            y2(1:lookback, :)=[];

            % Bollinger band strategy
            entryZscore=1;
            exitZscore=0;

            MA=movingAvg(yport, lookback);
            MSTD=movingStd(yport, lookback);
            zScore=(yport-MA)./MSTD;

            longsEntry=zScore < -entryZscore; % a long position means we should buy EWC
            longsExit=zScore > -exitZscore;

            shortsEntry=zScore > entryZscore;
            shortsExit=zScore < exitZscore;

            numUnitsLong=NaN(length(yport), 1);
            numUnitsShort=NaN(length(yport), 1);

            numUnitsLong(1)=0;
            numUnitsLong(longsEntry)=1; 
            numUnitsLong(longsExit)=0;
            numUnitsLong=fillMissingData(numUnitsLong); % fillMissingData can be downloaded from epchan.com/book2. It simply carry forward an existing position from previous day if today's positio is an indeterminate NaN.

            numUnitsShort(1)=0;
            numUnitsShort(shortsEntry)=-1; 
            numUnitsShort(shortsExit)=0;
            numUnitsShort=fillMissingData(numUnitsShort);

            numUnits=numUnitsLong+numUnitsShort;
            positions=repmat(numUnits, [1 size(y2, 2)]).*[-hedgeRatio ones(size(hedgeRatio))].*y2; % [hedgeRatio -ones(size(hedgeRatio))] is the shares allocation, [hedgeRatio -ones(size(hedgeRatio))].*y2 is the dollar capital allocation, while positions is the dollar capital in each ETF.
            pnl=sum(lag(positions, 1).*(y2-lag(y2, 1))./lag(y2, 1), 2); % daily P&L of the strategy
            ret=pnl./sum(abs(lag(positions, 1)), 2); % return is P&L divided by gross market value of portfolio
            ret(isnan(ret))=0;
            y_re = cumprod(1+ret)-1;
            h = figure_S53(y_re,tref,key_str);
            bac_re={y_re,key_str};            
        end
        function [H,re] = A38cis()
            key_str = 'S54-A38';
            sprintf('开始%s',key_str)
            index_pool = {'000016','000300';'000300','000905'};
            H=zeros(size(index_pool(:,1)));
            re = cell(size(index_pool(:,1)));
            for index_ID = 1:2
                sql_str = 'select tradeDate,closeIndex from yuqerdata.yq_index where symbol = "%s" order by tradeDate';
                sub_ticker1 = index_pool{index_ID,1};
                sub_ticker2 = index_pool{index_ID,2};
                x=fetchmysql(sprintf(sql_str,sub_ticker1),2);
                y=fetchmysql(sprintf(sql_str,sub_ticker2),2);
                [tref,ia,ib] = intersect(x(:,1),y(:,1));
                x = cell2mat(x(ia,2));
                y = cell2mat(y(ib,2));

                % Augment x with ones to  accomodate possible offset in the regression
                % between y vs x.

                x=[x ones(size(x))];

                %delta=0.0001; % delta=1 gives fastest change in beta, delta=0.000....1 allows no change (like traditional linear regression).
                delta = 1e-5;

                yhat=NaN(size(y)); % measurement prediction
                e=NaN(size(y)); % measurement prediction error
                Q=NaN(size(y)); % measurement prediction error variance

                % For clarity, we denote R(t|t) by P(t).
                % initialize R, P and beta.
                R=zeros(2);
                P=zeros(2);
                beta=NaN(2, size(x, 1));
                Vw=delta/(1-delta)*eye(2);
                Ve=0.001;


                % Initialize beta(:, 1) to zero
                beta(:, 1)=0;

                % Given initial beta and R (and P)
                for t=1:length(y)
                    if (t > 1)
                        beta(:, t)=beta(:, t-1); % state prediction. Equation 3.7
                        R=P+Vw; % state covariance prediction. Equation 3.8
                    end

                    yhat(t)=x(t, :)*beta(:, t); % measurement prediction. Equation 3.9

                    Q(t)=x(t, :)*R*x(t, :)'+Ve; % measurement variance prediction. Equation 3.10


                    % Observe y(t)
                    e(t)=y(t)-yhat(t); % measurement prediction error

                    K=R*x(t, :)'/Q(t); % Kalman gain

                    beta(:, t)=beta(:, t)+K*e(t); % State update. Equation 3.11
                    P=R-K*x(t, :)*R; % State covariance update. Euqation 3.12

                end
                y2=[x(:, 1) y];
                %adair modified
                longsEntry=e > -sqrt(Q); % a long position means we should buy EWC
                longsExit=e < -sqrt(Q);

                shortsEntry=e < sqrt(Q);
                shortsExit=e > sqrt(Q);

                numUnitsLong=NaN(length(y2), 1);
                numUnitsShort=NaN(length(y2), 1);

                numUnitsLong(1)=0;
                numUnitsLong(longsEntry)=1; 
                numUnitsLong(longsExit)=0;
                numUnitsLong=fillMissingData(numUnitsLong); % fillMissingData can be downloaded from epchan.com/book2. It simply carry forward an existing position from previous day if today's positio is an indeterminate NaN.

                numUnitsShort(1)=0;
                numUnitsShort(shortsEntry)=-1; 
                numUnitsShort(shortsExit)=0;
                numUnitsShort=fillMissingData(numUnitsShort);

                numUnits=numUnitsLong+numUnitsShort;
                positions=repmat(numUnits, [1 size(y2, 2)]).*[-beta(1, :)' ones(size(beta(1, :)'))].*y2; % [hedgeRatio -ones(size(hedgeRatio))] is the shares allocation, [hedgeRatio -ones(size(hedgeRatio))].*y2 is the dollar capital allocation, while positions is the dollar capital in each ETF.
                pnl=sum(lag(positions, 1).*(y2-lag(y2, 1))./lag(y2, 1), 2); % daily P&L of the strategy

                ret=pnl./sum(abs(lag(positions, 1)), 2); % return is P&L divided by gross market value of portfolio
                ret(isnan(ret))=0;
                y_re = cumprod(1+ret)-1;
                h = figure_S53(y_re,tref,[]);
                t_str = sprintf('%s-%s-%s',key_str,sub_ticker1,sub_ticker2);
                title(t_str)
                H(index_ID) = h;
                re{index_ID} = {y_re,t_str}';
            end
            re = [re{:}];
            
        end
        function [h,bac_re] = A38_EWAEWC()
            key_str = 'S54-A38EWAEWC';
            sprintf('开始%s',key_str)
            sql_str = 'select tradeDate,closePrice from polygon.usastock_day where ticker = "%s" order by tradeDate';
            x=fetchmysql(sprintf(sql_str,'EWA'),2);
            y=fetchmysql(sprintf(sql_str,'EWC'),2);
            [tref,ia,ib] = intersect(x(:,1),y(:,1));
            x = cell2mat(x(ia,2));
            y = cell2mat(y(ib,2));

            % Augment x with ones to  accomodate possible offset in the regression
            % between y vs x.

            x=[x ones(size(x))];

            delta=0.0001; % delta=1 gives fastest change in beta, delta=0.000....1 allows no change (like traditional linear regression).

            yhat=NaN(size(y)); % measurement prediction
            e=NaN(size(y)); % measurement prediction error
            Q=NaN(size(y)); % measurement prediction error variance

            % For clarity, we denote R(t|t) by P(t).
            % initialize R, P and beta.
            R=zeros(2);
            P=zeros(2);
            beta=NaN(2, size(x, 1));
            Vw=delta/(1-delta)*eye(2);
            Ve=0.001;
            % Initialize beta(:, 1) to zero
            beta(:, 1)=0;

            % Given initial beta and R (and P)
            for t=1:length(y)
                if (t > 1)
                    beta(:, t)=beta(:, t-1); % state prediction. Equation 3.7
                    R=P+Vw; % state covariance prediction. Equation 3.8
                end

                yhat(t)=x(t, :)*beta(:, t); % measurement prediction. Equation 3.9

                Q(t)=x(t, :)*R*x(t, :)'+Ve; % measurement variance prediction. Equation 3.10


                % Observe y(t)
                e(t)=y(t)-yhat(t); % measurement prediction error

                K=R*x(t, :)'/Q(t); % Kalman gain

                beta(:, t)=beta(:, t)+K*e(t); % State update. Equation 3.11
                P=R-K*x(t, :)*R; % State covariance update. Euqation 3.12

            end

            y2=[x(:, 1) y];
            longsEntry=e < -sqrt(Q); % a long position means we should buy EWC
            longsExit=e > -sqrt(Q);

            shortsEntry=e > sqrt(Q);
            shortsExit=e < sqrt(Q);

            numUnitsLong=NaN(length(y2), 1);
            numUnitsShort=NaN(length(y2), 1);

            numUnitsLong(1)=0;
            numUnitsLong(longsEntry)=1; 
            numUnitsLong(longsExit)=0;
            numUnitsLong=fillMissingData(numUnitsLong); % fillMissingData can be downloaded from epchan.com/book2. It simply carry forward an existing position from previous day if today's positio is an indeterminate NaN.

            numUnitsShort(1)=0;
            numUnitsShort(shortsEntry)=-1; 
            numUnitsShort(shortsExit)=0;
            numUnitsShort=fillMissingData(numUnitsShort);

            numUnits=numUnitsLong+numUnitsShort;
            positions=repmat(numUnits, [1 size(y2, 2)]).*[-beta(1, :)' ones(size(beta(1, :)'))].*y2; % [hedgeRatio -ones(size(hedgeRatio))] is the shares allocation, [hedgeRatio -ones(size(hedgeRatio))].*y2 is the dollar capital allocation, while positions is the dollar capital in each ETF.
            pnl=sum(lag(positions, 1).*(y2-lag(y2, 1))./lag(y2, 1), 2); % daily P&L of the strategy
            ret=pnl./sum(abs(lag(positions, 1)), 2); % return is P&L divided by gross market value of portfolio
            ret(isnan(ret))=0;

            y_re = cumprod(1+ret)-1;
            %setfigure
            t_str = 'A38-EWAEWC';
            h = figure_S53(y_re,tref,t_str);
            %title(t_str);
            bac_re = {y_re,t_str};            
        end
        function [H,re] = A38_forex()
            key_str = 'S54-A38forex';
            sprintf('开始%s',key_str)
            sql_str = 'select tradeDate,closePrice from aksharedata.currency_hist where ticker = "%s" and tradeDate>="2011-01-01" order by tradeDate';
            %usdtwd usdcnh
            symbol_pool = {'usdtwd','usdkrw';'usdtwd','usdcnh';'usdtry','usdzar'};
            delta_pool = [1e-4,1e-7];
            H = zeros(2,1);
            re = cell(2,1);
            for pool_id = 1:2
                sym1 = symbol_pool{pool_id,1};
                sym2 = symbol_pool{pool_id,2};

                x=fetchmysql(sprintf(sql_str,sym1),2);
                y=fetchmysql(sprintf(sql_str,sym2),2);

                [tref,ia,ib] = intersect(x(:,1),y(:,1));
                x = cell2mat(x(ia,2));
                y = cell2mat(y(ib,2));

                % Augment x with ones to  accomodate possible offset in the regression
                % between y vs x.

                x=[x ones(size(x))];

                %delta=0.0001; % delta=1 gives fastest change in beta, delta=0.000....1 allows no change (like traditional linear regression).
                %delta = 1e-7;
                %delta = 1e-4;
                delta = delta_pool(pool_id);
                yhat=NaN(size(y)); % measurement prediction
                e=NaN(size(y)); % measurement prediction error
                Q=NaN(size(y)); % measurement prediction error variance

                % For clarity, we denote R(t|t) by P(t).
                % initialize R, P and beta.
                R=zeros(2);
                P=zeros(2);
                beta=NaN(2, size(x, 1));
                Vw=delta/(1-delta)*eye(2);
                Ve=0.001;


                % Initialize beta(:, 1) to zero
                beta(:, 1)=0;

                % Given initial beta and R (and P)
                for t=1:length(y)
                    if (t > 1)
                        beta(:, t)=beta(:, t-1); % state prediction. Equation 3.7
                        R=P+Vw; % state covariance prediction. Equation 3.8
                    end

                    yhat(t)=x(t, :)*beta(:, t); % measurement prediction. Equation 3.9

                    Q(t)=x(t, :)*R*x(t, :)'+Ve; % measurement variance prediction. Equation 3.10


                    % Observe y(t)
                    e(t)=y(t)-yhat(t); % measurement prediction error

                    K=R*x(t, :)'/Q(t); % Kalman gain

                    beta(:, t)=beta(:, t)+K*e(t); % State update. Equation 3.11
                    P=R-K*x(t, :)*R; % State covariance update. Euqation 3.12

                end

                %{
                plot(beta(1, :)');
                figure;
                plot(beta(2, :)');
                figure;
                plot(e(3:end), 'r');
                hold on;
                plot(sqrt(Q(3:end)));
                %}
                y2=[x(:, 1) y];
                longsEntry=e < -sqrt(Q); % a long position means we should buy EWC
                longsExit=e > -sqrt(Q);

                shortsEntry=e > sqrt(Q);
                shortsExit=e < sqrt(Q);

                numUnitsLong=NaN(length(y2), 1);
                numUnitsShort=NaN(length(y2), 1);

                numUnitsLong(1)=0;
                numUnitsLong(longsEntry)=1; 
                numUnitsLong(longsExit)=0;
                numUnitsLong=fillMissingData(numUnitsLong); % fillMissingData can be downloaded from epchan.com/book2. It simply carry forward an existing position from previous day if today's positio is an indeterminate NaN.

                numUnitsShort(1)=0;
                numUnitsShort(shortsEntry)=-1; 
                numUnitsShort(shortsExit)=0;
                numUnitsShort=fillMissingData(numUnitsShort);

                numUnits=numUnitsLong+numUnitsShort;
                positions=repmat(numUnits, [1 size(y2, 2)]).*[-beta(1, :)' ones(size(beta(1, :)'))].*y2; % [hedgeRatio -ones(size(hedgeRatio))] is the shares allocation, [hedgeRatio -ones(size(hedgeRatio))].*y2 is the dollar capital allocation, while positions is the dollar capital in each ETF.
                pnl=sum(lag(positions, 1).*(y2-lag(y2, 1))./lag(y2, 1), 2); % daily P&L of the strategy
                ret=pnl./sum(abs(lag(positions, 1)), 2); % return is P&L divided by gross market value of portfolio
                ret(isnan(ret))=0;
                t_str = sprintf('A38-%s vs %s',sym1,sym2);
                y_re = cumprod(1+ret)-1;
                %setfigure
                H(pool_id) = figure_S53(y_re,tref,t_str);
                re{pool_id}={y_re,t_str}';                
            end
            re = [re{:}];
        end
        function [h,re] = A41_bog()
            key_str = 'A41-bog';
            sprintf('开始%s',key_str)
            topN=10; % Max number of positions
            entryZscore=1;
            lookback=20; % for MA
            %load('../Data/inputDataOHLCDaily_20120424', 'stocks', 'tday', 'op', 'hi', 'lo', 'cl');
            %load('inputDataOHLCDaily_stocks_20120424', 'stocks', 'tday', 'op', 'hi', 'lo', 'cl');
            key_str = 'A_41';
            sql_str = 'select * from S54.A41data order by tradeDate';
            X = fetchmysql(sql_str,2);
            %tref = cellfun(@(x) x(:,1)',X,'UniformOutput',false);
            tref = unique(X(:,2));
            tday = cellfun(@(x) str2double(strjoin(strsplit(x,'-'))),tref);
            stocks = unique(X(:,1));
            T = length(stocks);
            x=cell(T,1);

            parfor i = 1:T
                ind = strcmp(X(:,1),stocks(i));
                temp = X(ind,2:end);
                if ~isempty(temp)
                    [~,ia,ib] = intersect(tref,temp(:,1));
                    sub_x = nan(length(tref),size(temp,2)-1);
                    sub_x(ia,:) = cell2mat(temp(ib,2:end));
                else
                    keyboard
                end
                x{i} = sub_x;
            end

            temp = cell(4,1);
            for i = 1:4
                sub_temp = cellfun(@(x) x(:,i),x,'UniformOutput',false);
                temp{i} = [sub_temp{:}];
            end

            op=temp{1};
            hi=temp{2};
            lo=temp{3};
            cl=temp{4};

            stdretC2C90d=backshift(1, smartMovingStd(calculateReturns(cl, 1), 90));
            buyPrice=backshift(1, lo).*(1-entryZscore*stdretC2C90d);

            retGap=(op-backshift(1, lo))./backshift(1, lo);

            pnl=zeros(length(tday), 1);

            positionTable=zeros(size(cl));

            ma=backshift(1, smartMovingAvg(cl, lookback));

            for t=2:size(cl, 1)
                hasData=find(isfinite(retGap(t, :)) & op(t, :) < buyPrice(t, :) & op(t, :) > ma(t, :));

                [foo idxSort]=sort(retGap(t, hasData), 'ascend');
                positionTable(t, hasData(idxSort(1:min(topN, length(idxSort)))))=1;
            end

            retO2C=(cl-op)./op;
            pnl=smartsum(positionTable.*(retO2C), 2);
            ret=pnl/topN; 
            ret(isnan(ret))=0;

            %fprintf(1, '%i - %i\n', tday(1), tday(end));
            %fprintf(1, 'APR=%10.4f\n', prod(1+ret).^(252/length(ret))-1);

            %fprintf(1, 'Sharpe=%4.2f\n', mean(ret)*sqrt(252)/std(ret));
            % APR=8.7%, Sharpe=1.5

            %cumret=cumprod(1+ret)-1; % compounded ROE

            %plot(cumret);
            y_re = cumprod(1+ret)-1;
            %setfigure
            t_str = sprintf('%s-%s',key_str,'bog');
            h = figure_S53(y_re,tref,t_str);
            re = {y_re,t_str};

        end
        function [h,re] = A41_csi()
            key_str = 'A41-CSI';
            sprintf('开始%s',key_str)
            topN=10; % Max number of positions
            entryZscore=1;
            lookback=20; % for MA
            index_pool = {'000300','000905','000001'};
            h = zeros(size(index_pool));
            re = cell(size(h));
            for index_id = 1:3
                index_sel = index_pool{index_id};
                stocks = yq_methods.get_index_pool(index_sel,'2005-01-01');

                temp = get_pool_data(stocks,'2005-01-01','csi');
                cl=temp.cl;
                hi=temp.hi;
                lo=temp.lo;
                op=temp.op;
                tday=temp.tday;
                tref = temp.tref;


                stdretC2C90d=backshift(1, smartMovingStd(calculateReturns(cl, 1), 90));
                buyPrice=backshift(1, lo).*(1-entryZscore*stdretC2C90d);

                retGap=(op-backshift(1, lo))./backshift(1, lo);

                pnl=zeros(length(tday), 1);

                positionTable=zeros(size(cl));

                ma=backshift(1, smartMovingAvg(cl, lookback));

                for t=2:size(cl, 1)
                    hasData=find(isfinite(retGap(t, :)) & op(t, :) < buyPrice(t, :) & op(t, :) > ma(t, :));

                    [foo idxSort]=sort(retGap(t, hasData), 'ascend');
                    positionTable(t, hasData(idxSort(1:min(topN, length(idxSort)))))=1;
                end

                retO2C=(cl-op)./op;


                pnl=smartsum(positionTable.*(retO2C), 2);
                ret=pnl/topN; 
                ret(isnan(ret))=0;

                %fprintf(1, '%i - %i\n', tday(1), tday(end));
                %fprintf(1, 'APR=%10.4f\n', prod(1+ret).^(252/length(ret))-1);

                %fprintf(1, 'Sharpe=%4.2f\n', mean(ret)*sqrt(252)/std(ret));
                % APR=8.7%, Sharpe=1.5

                %cumret=cumprod(1+ret)-1; % compounded ROE

                %plot(cumret);
                y_re = cumprod(1+ret)-1;
                %setfigure
                t_str = sprintf('%s-%s',key_str,index_sel);
                h(index_id) = figure_S53(y_re,tref,t_str);
                re{index_id} = {y_re,t_str}';
            end
            re = [re{:}];
            
        end
        
        function [h,re] = A43_indexArb()
            key_str = 'S54-A43';
            sprintf('开始%s',key_str)
            etf0=[];
            %stks = get_pool_data(stks0.stocks);
            sql_str = 'select * from S54.A41data order by tradeDate';
            X = fetchmysql(sql_str,2);
            %tref = cellfun(@(x) x(:,1)',X,'UniformOutput',false);
            tref = unique(X(:,2));
            tday = cellfun(@(x) str2double(strjoin(strsplit(x,'-'),'')),tref);
            stocks = unique(X(:,1));
            T = length(stocks);
            x=cell(T,1);
            parfor i = 1:T
                ind = strcmp(X(:,1),stocks(i));
                temp = X(ind,2:end);
                if ~isempty(temp)
                    [~,ia,ib] = intersect(tref,temp(:,1));
                    sub_x = nan(length(tref),size(temp,2)-1);
                    sub_x(ia,:) = cell2mat(temp(ib,2:end));
                else
                    keyboard
                end
                x{i} = sub_x;
            end
            temp = cell(4,1);
            for i = 1:4
                sub_temp = cellfun(@(x) x(:,i),x,'UniformOutput',false);
                temp{i} = [sub_temp{:}];
            end
            stks=[];
            stks.op=temp{1};
            stks.hi=temp{2};
            stks.lo=temp{3};
            stks.cl=temp{4};
            stks.stocks= stocks;
            stks.tday=tday;
            stks.tref = tref;
            etf0.syms = {'SPY'};
            etf = get_pool_data(etf0.syms);
            etf.syms = etf.stocks;
            %}

            % Ensure data have same dates
            [tday,idx1,idx2]=intersect(stks.tday, etf.tday);

            tref = stks.tref(idx1);

            stks.cl=stks.cl(idx1, :);
            etf.cl=etf.cl(idx2, :);

            % Use SPY
            idxS=find(strcmp('SPY', etf.syms));
            etf.cl=etf.cl(:, idxS);

            trainDataIdx=find(tday>=20070101 & tday<=20071231);
            testDataIdx=find(tday > 20071231);
            if ~isempty(tref)
                tref = tref(tday>20071231);
            end

            isCoint=false(size(stks.stocks));
            for s=1:length(stks.stocks)
                % Combine the two time series into a matrix y2 for input into Johansen test
                y2=[stks.cl(trainDataIdx, s), etf.cl(trainDataIdx)];
                badData=any(isnan(y2), 2);
                y2(badData, :)=[]; % remove any missing data

                if (size(y2, 1) > 250)
                    results=johansen(y2, 0, 1); % johansen test with non-zero offset but zero drift, and with the lag k=1.
                    if (results.lr1(1) > results.cvt(1, 1))
                        isCoint(s)=true;
                    end
                end    
            end

            length(find(isCoint))
            % 98: there are 98 stocks that are cointegrating with SPY

            % Form a long-only portfolio with all stocks that cointegrate with SPY, with equal
            % capital allocation
            yN=stks.cl(trainDataIdx, isCoint);
            logMktVal_long=sum(log(yN), 2); % The net market value of the long-only portfolio is same as the "spread"

            % Confirm that the portfolio cointegrates with SPY
            ytest=[logMktVal_long, log(etf.cl(trainDataIdx))]; 
            results=johansen(ytest, 0, 1); % johansen test with non-zero offset but zero drift, and with the lag k=1.
            prt(results);           
            
            % Apply linear mean-reversion model on test set
            yNplus=[stks.cl(testDataIdx, isCoint), etf.cl(testDataIdx)]; % Array of stock and ETF prices
            weights=[repmat(results.evec(1, 1), size(stks.cl(testDataIdx, isCoint))), ...
                   repmat(results.evec(2, 1), size(etf.cl(testDataIdx)))]; % Array of log market value of stocks and ETF's

            logMktVal=smartsum(weights.*log(yNplus), 2); % Log market value of long-short portfolio

            lookback=5;
            numUnits=-(logMktVal-movingAvg(logMktVal, lookback))./movingStd(logMktVal, lookback); % capital invested in portfolio in dollars.  movingAvg and movingStd are functions from epchan.com/book2
            positions=repmat(numUnits, [1 size(weights, 2)]).*weights; % positions is the dollar capital in each stock or ETF.
            pnl=smartsum(lag(positions, 1).*(log(yNplus)-lag(log(yNplus), 1)), 2); % daily P&L of the strategy
            ret=pnl./smartsum(abs(lag(positions, 1)), 2); % return is P&L divided by gross market value of portfolio
            ret(isnan(ret))=0;

            %figure;
            %plot(cumprod(1+ret)-1); % Cumulative compounded return
            y_re = cumprod(1+ret)-1;
            %setfigure
            t_str = sprintf('%s-indexArb',key_str);
            h = figure_S53(y_re,tref,t_str);
            re = {y_re,t_str};
        end
        
        function [h,re] = A44_andrewlo()
            key_str = 'S54-A44andrewlo';
            sprintf('开始%s',key_str)
            h = zeros(2,1);
            re = cell(size(h));
            stks = get_A41data();
            %stks = get_pool_data_update(stocks);
            tday = stks.tday;
            cl=stks.cl;
            op=stks.op;
            tref = stks.tref;
            ind = tday>=20070103;
            %ind = tday>=20070103 & tday<=20111230;
            tday=tday(ind);
            cl=cl(ind, :);
            op=op(ind, :);
            if ~isempty(tref)
            tref=tref(ind);
            end
            % cl is a TxN array of closing prices, where T is the number of trading
            % days, and N is the number of stocks in the S&P 500
            ret=(cl-lag(cl, 1))./lag(cl, 1); % daily returns

            marketRet=smartmean(ret, 2); % equal weighted market index return

            weights=-(ret-repmat(marketRet, [1 size(ret, 2)]));
            weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

            dailyret=smartsum(backshift(1, weights).*ret, 2); % Capital is always one

            dailyret(isnan(dailyret))=0;

            %plot(cumprod(1+dailyret)-1); % Cumulative compounded return
            y_re = cumprod(1+dailyret)-1;
            %setfigure
            t_str = sprintf('%s-closeprice-trade',key_str);
            h(1) = figure_S53(y_re,tref,t_str);
            re{1} = {y_re,t_str}';

            %fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+dailyret).^(252/length(dailyret))-1, sqrt(252)*mean(dailyret)/std(dailyret));
            % APR=13.7%, Sharpe=1.3

            % daily pnl with transaction costs deducted
            % onewaytcost=0.0005; % assume 5 basis points
            % 
            % dailyretMinustcost=dailyret - ...
            %     smartsum(abs(weights./cl-backshift(1, weights)./backshift(1, cl)).*backshift(1, cl), 2).*onewaytcost./smartsum(abs(weights), 2); % transaction costs are only incurred when the weights change
            % 
            % annavgretMinustcost=252*smartmean(dailyretMinustcost, 1)*100
            % 
            % sharpeMinustcost=sqrt(252)*smartmean(dailyretMinustcost, 1)/smartstd(dailyretMinustcost, 1) 
            % 
            % % switch to use open prices
            % 
            ret=(op-backshift(1, cl))./backshift(1, cl); % daily returns

            marketRet=smartmean(ret, 2); % equal weighted market index return

            weights=-(ret-repmat(marketRet, [1 size(ret, 2)])); % weight of a stock is proportional to the negative distance to the market index.
            weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

            dailyret=smartsum(weights.*(cl-op)./op, 2)./smartsum(abs(weights), 2);
            dailyret(isnan(dailyret))=0;

            %plot(cumprod(1+dailyret)-1); % Cumulative compounded return
            y_re = cumprod(1+dailyret)-1;
            %setfigure
            t_str = sprintf('%s-openprice-trade',key_str);
            h(2) = figure_S53(y_re,tref,t_str);
            re{2} = {y_re,t_str}';
            
            re = [re{:}];
            
        end
        function [h,re] = A44_csi()
            key_str = 'S54-A44csi';
            sprintf('开始%s',key_str)
            index_pool = {'000300','000905','000001'};
            fee = 3/1000;
            h=zeros(3,1);
            re = cell(size(h));
            for index_id = 1:3
                index_sel = index_pool{index_id};
                stocks = yq_methods.get_index_pool(index_sel,'2005-01-01');
                stks = get_pool_data_update(stocks,'2005-01-01','csi');
                %load tempstks.mat
                tday = stks.tday;
                cl=stks.cl;
                op=stks.op;
                tref = stks.tref;
                %idxStart=find(tday==20070103);
                %idxEnd=find(tday==20111230);
                ind = tday>=20100101;
                %ind = tday>=20070103 & tday<=20111230;
                tday=tday(ind);
                cl=cl(ind, :);
                op=op(ind, :);
                tref=tref(ind);
                % cl is a TxN array of closing prices, where T is the number of trading
                % days, and N is the number of stocks in the S&P 500
                ret=(cl-lag(cl, 1))./lag(cl, 1); % daily returns

                marketRet=smartmean(ret, 2); % equal weighted market index return

                weights=-(ret-repmat(marketRet, [1 size(ret, 2)]));
                weights(weights<0)=0;
                weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

                dailyret=smartsum(backshift(1, weights).*ret, 2); % Capital is always one

                dailyret(isnan(dailyret))=0;

%                 %plot(cumprod(1+dailyret)-1); % Cumulative compounded return
%                 if ~isempty(tref)
%                     y_re = cumprod(1+dailyret-fee)-1;
%                     %setfigure
%                     h = figure_S53(y_re,tref,[]);
%                     title(sprintf('%s-%s-part1',key_str,index_sel))
%                 else
%                     figure;
%                     plot(cumprod(1+dailyret)-1); % Cumulative compounded return
%                 end

                %fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+dailyret).^(252/length(dailyret))-1, sqrt(252)*mean(dailyret)/std(dailyret));
                % APR=13.7%, Sharpe=1.3

                % daily pnl with transaction costs deducted
                % onewaytcost=0.0005; % assume 5 basis points
                % 
                % dailyretMinustcost=dailyret - ...
                %     smartsum(abs(weights./cl-backshift(1, weights)./backshift(1, cl)).*backshift(1, cl), 2).*onewaytcost./smartsum(abs(weights), 2); % transaction costs are only incurred when the weights change
                % 
                % annavgretMinustcost=252*smartmean(dailyretMinustcost, 1)*100
                % 
                % sharpeMinustcost=sqrt(252)*smartmean(dailyretMinustcost, 1)/smartstd(dailyretMinustcost, 1) 
                % 
                % % switch to use open prices
                % 
                ret=(op-backshift(1, cl))./backshift(1, cl); % daily returns

                marketRet=smartmean(ret, 2); % equal weighted market index return

                weights=-(ret-repmat(marketRet, [1 size(ret, 2)])); % weight of a stock is proportional to the negative distance to the market index.
                weights(weights<0)=0;
                weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

                dailyret=smartsum(weights.*(cl-op)./op, 2)./smartsum(abs(weights), 2);
                dailyret(isnan(dailyret))=0;

                %plot(cumprod(1+dailyret)-1); % Cumulative compounded return

                y_re = cumprod(1+dailyret-fee)-1;
                %setfigure
                t_str = sprintf('%s-%s-openprice-trade',key_str,index_sel);
                h(index_id) = figure_S53(y_re,tref,t_str);
                re{index_id} = {y_re,t_str}';
            end
            re = [re{:}]; 
        end
        function [h,re] = A51_audcad()
            key_str = 'S54-A51'; 
            sprintf('开始%s',key_str)
            sql_str = ['select date(tradeDate),closePrice from polygon_fx_minute.%s ',...
                'where tradeDate>="2009-01-01" and 100*hour(tradeDate)+minute(tradeDate)=359 order by tradeDate'];
            usdcad = fetchmysql(sprintf(sql_str,'USDCAD'),2);
            audusd = fetchmysql(sprintf(sql_str,'AUDUSD'),2);
            [tref,ia,ib] = intersect(usdcad(:,1),audusd(:,1));
            usdcad_cl=cell2mat(usdcad(ia,2));
            audusd_cl = cell2mat(audusd(ib,2));
            tday = cellfun(@(x) str2double(strjoin(strsplit(x,'-'),'')),tref);

            % Need to invert currency pair so that each unit has same capital in local
            % currency
            cad=1./usdcad_cl;
            aud=audusd_cl;

            y=[ aud,cad   ];
            trainlen=250;
            lookback=20;
            hedgeRatio=NaN(size(y));
            numUnits=NaN(size(y, 1), 1);

            for t=trainlen+1:size(y, 1)
                res=johansen(y(t-trainlen:t-1, :), 0, 1);
                hedgeRatio(t, :)=res.evec(:, 1)';

                % yport is the market value of a unit portfolio of AUDUSD and CADUSD expressed in US$.
                yport=sum(y(t-lookback+1:t, :).*repmat(hedgeRatio(t, :), [lookback 1]), 2);
                ma=mean(yport);
                mstd=std(yport);
                zScore=(yport(end)-ma)/mstd;

                % numUnits are number of units of unit portfolio of AUDUSD and CADUSD
                numUnits(t)=-(yport(end)-ma)/mstd;

            end
            tref = tref(trainlen+1:end);

            % positions are market values of AUDUSD and CADUSD in portfolio expressed
            % in US$.
            positions=repmat(numUnits, [1 size(y, 2)]).*hedgeRatio.*y; 

            % daily P&L of portfolio in US$.
            pnl=sum(lag(positions, 1).*(y-lag(y, 1))./lag(y, 1), 2); 
            ret=pnl./sum(abs(lag(positions, 1)), 2); 
            ret(isnan(ret))=0;

            %plot(cumprod(1+ret(trainlen+1:end))-1); % Cumulative compounded return

            y_re = cumprod(1+ret(trainlen+1:end))-1;
            %setfigure
            t_str = sprintf('%s-AUDCAD',key_str);
            h = figure_S53(y_re,tref,t_str);
            re = {y_re,t_str};
        end
        function [h,re] = A57()
            key_str = 'S54-A57';
            sprintf('开始%s',key_str)
            % load('inputDataDaily_VX_20120507', 'tday', 'contracts', 'cl');
            %load('inputDataDaily_CL_20120813', 'tday', 'contracts', 'cl');
            ticker_pool = {'CL','TU','HG','PA'};
            h = zeros(size(ticker_pool));
            re = cell(size(h));
            for ticker_ID = 1:length(ticker_pool)
                %ticker = 'TU';
                ticker = ticker_pool{ticker_ID};
                sql_str = ['select year(tradeDate)*10000+month(tradeDate)*100+day(tradeDate) as t, ',...
                    'ticker,settlePrice from S50.us_futures where mark4 = "%s" and tradeDate>="2006-01-01" order by tradeDate'];
                X = fetchmysql(sprintf(sql_str,ticker),3);
                % T=[1:length(spot)]';
                % isBadData=~isfinite(spot);
                % spot(isBadData)=[];
                % T(isBadData)=[];
                % res=ols(log(spot), [T ones(size(T, 1), 1)]);
                % 
                % fprintf(1, 'Average annualized spot return=%f\n', 252*smartmean(res.beta(1)));
                %f1 = @(x) [x(end-3:end),x(1:end-4)];
                a=table2cell(X(:,2));
                b=cellfun(@(x)  [x(end-3:end),x(1:end-4)],a,'UniformOutput',false);
                X.ticker=b;
                %X.ticker = varfun(@(x) [x(end-3:end),x(1:end-4)],X.ticker);
                X = unstack(X,'settlePrice','ticker');

                cl = table2array(X(:,2:end));
                tday = X.t;
                contracts = X.Properties.VariableNames(2:end);
                % Fitting gamma to forward curve
                gamma=NaN(size(tday));
                for t=1:length(tday)

                    FT=cl(t, :)';
                    idx=find(isfinite(FT));
                    idxDiff=fwdshift(1, idx)-idx; % ensure consecutive months futures
                    if (length(idx) >= 5 && all(idxDiff(1:4)==1))
                        FT=FT(idx(1:5)); % only uses the nearest 5 contracts
                        T=[1:length(FT)]';
                %         scatter(T, log(FT));
                        res=ols(log(FT), [T ones(size(T, 1), 1)]);
                        gamma(t)=-12*res.beta(1);
                    end
                end
                gamma=fillMissingData(gamma);

                % plot(gamma);


                %print -r300 -djpeg fig5_4
                % hold on;

                % fprintf(1, 'Average annualized roll return=%f\n', smartmean(gamma));
                isGoodData=find(isfinite(gamma));
                results=adf(gamma(isGoodData), 0, 1);
                prt(results);

                gammalag=lag(gamma(isGoodData), 1);  
                deltaGamma=gamma(isGoodData)-gammalag;
                deltaGamma(1)=[]; 
                gammalag(1)=[];
                regress_results=ols(deltaGamma, [gammalag ones(size(gammalag))]);
                halflife=-log(2)/regress_results.beta(1);

                %fprintf(1, 'halflife=%f days\n', halflife);
                % halflife=36.394034 days


                lookback=round(halflife);
                ma=movingAvg(gamma, lookback);
                mstd=movingStd(gamma, lookback);
                zScore=(gamma-ma)./mstd;

                % linear mean reversion strategy
                isExpireDate=false(size(cl));
                positions=zeros(size(cl));

                isExpireDate=isfinite(cl) & ~isfinite(fwdshift(1, cl));

                holddays=3*21;
                numDaysStart=holddays+10;
                numDaysEnd=10;
                spreadMonth=12; % No. months between far and near contracts.
                for c=1:length(contracts)-spreadMonth
                    expireIdx=find(isExpireDate(:, c));
                    expireIdx=expireIdx(end); % There may be some missing data earlier on
                    if (c==1)
                        startIdx=max(1, expireIdx-numDaysStart);
                        endIdx=expireIdx-numDaysEnd;
                    else % ensure next front month contract doesn't start until current one ends
                        myStartIdx=endIdx+1;
                        myEndIdx=expireIdx-numDaysEnd;
                        if (myEndIdx-myStartIdx >= holddays)
                            startIdx=myStartIdx;
                            endIdx=myEndIdx;
                        else
                            startIdx=NaN;
                        end
                    end

                    if (~isempty(expireIdx) && endIdx > startIdx)
                        positions(startIdx:endIdx, c)=-1; % Presume we long spread (long back contract, short front contract)
                        positions(startIdx:endIdx, c+spreadMonth)=1;
                    end
                end



                positions(isnan(zScore), :)=0;
                positions(zScore > 0, :)=-positions(zScore > 0, :);
                % positions(zScore > 1, :)=-positions(zScore > 1, :);

                ret=smartsum(lag(positions).*(cl-lag(cl, 1))./lag(cl, 1), 2)/2;
                ret(isnan(ret))=0;

                idx=find(tday==20080102);
                % idx=1;

                cumret=cumprod(1+ret(idx:end))-1;
                y_re = cumret; 
                tref = cellfun(@num2str,num2cell(tday(idx:end)),'UniformOutput',false);
                tref = datenum(tref,'yyyymmdd');
                tref = cellstr(datestr(tref,'yyyy-mm-dd'));
                %setfigure
                t_str = sprintf('%s-%s',key_str,ticker);
                h(ticker_ID) = figure_S53(y_re,tref,t_str);
                re{ticker_ID} = {y_re,t_str}';
            end
            re = [re{:}];
            
        end
        function [h,re] = A63()
            key_str = 'S54-A63';
            sprintf('开始%s',key_str)
            temp = get_pool_data({'USO','XLE'});
            uso=temp.cl(:,1);
            xle = temp.cl(:,2);
            tday_ETF = temp.tday;

            ticker = 'PA';
            sql_str = ['select year(tradeDate)*10000+month(tradeDate)*100+day(tradeDate) as t, ',...
                'ticker,settlePrice from S50.us_futures where mark4 = "%s" and tradeDate>="2000-11-20"  ',...
                'and abs(year(tradeDate)-mark2)<=3 order by tradeDate'];
            X = fetchmysql(sprintf(sql_str,ticker),3);
            a=table2cell(X(:,2));
            b=cellfun(@(x)  [x(end-3:end),x(1:end-4)],a,'UniformOutput',false);
            X.ticker=b;
            %X.ticker = varfun(@(x) [x(end-3:end),x(1:end-4)],X.ticker);
            X = unstack(X,'settlePrice','ticker');
            cl = table2array(X(:,2:end));
            tday = X.t;
            contracts = X.Properties.VariableNames(2:end);
            ratioMatrix=(fwdshift(1, cl')./cl')'; % back/front

            ratio=NaN(size(ratioMatrix, 1), 1);
            isExpireDate=false(size(ratio));

            isExpireDate=isfinite(cl) & ~isfinite(fwdshift(1, cl));

            % Define front month as 40 days to 10 days before expiration
            numDaysStart=40;
            numDaysEnd=10;

            for c=1:length(contracts)-1
                expireIdx=find(isExpireDate(:, c),1,'last');
                if (c==1)
                    startIdx=expireIdx-numDaysStart;
                    endIdx=expireIdx-numDaysEnd;
                else % ensure next front month contract doesn't start until current one ends
                    startIdx=max(endIdx+1, expireIdx-numDaysStart);
                    endIdx=expireIdx-numDaysEnd;
                end

                if (~isempty(expireIdx)) && startIdx>0
                    ratio(startIdx:endIdx)=ratioMatrix(startIdx:endIdx, c);
                end
            end

            [tday,idxA,idxB]=intersect(tday_ETF, tday);
            uso=uso(idxA);
            xle=xle(idxA);
            ratio=ratio(idxB);

            positions=zeros(length(tday), 2);

            % Contango, negative roll return, buy spot, short future
            contango=find(ratio > 1);
            positions(contango, :)=repmat([-1 1], [length(contango) 1]);
            % Backwardation, positive roll return, short spot, long future
            backwardation=find(ratio < 1);
            positions(backwardation, :)=repmat([1 -1], [length(backwardation) 1]);

            ret=smartsum(lag(positions, 1).*([uso xle]-lag([uso xle], 1))./lag([uso xle], 1), 2)/2;

            ret(isnan(ret))=0;

            cumret=cumprod(1+ret)-1;
            y_re = cumret; 
            tref = cellfun(@num2str,num2cell(tday),'UniformOutput',false);
            tref = datenum(tref,'yyyymmdd');
            tref = cellstr(datestr(tref,'yyyy-mm-dd'));
            %setfigure
            t_str = sprintf('%s-%s',key_str,ticker);
            h = figure_S53(y_re,tref,t_str);
            re = {y_re,t_str};

        end
        function [h,re] = A71()
            key_str = 'S54-A71';
            sprintf('开始%s',key_str)
            %update data update_A71data
            update_A71data();
            entryZscore=0.1;
            symbol = 'FSTX';
            data = fetchmysql('select * from S54.A71data',2);
            %[~,~,data]= xlsread('a71-2.csv');
            %data = data(2:end,:);
            tref_num = datenum(data(:,1));
            [tref_num,ia] = sort(tref_num);
            tref = cellstr(datestr(tref_num,'yyyy-mm-dd'));

            %id=tref_num<datenum(2012,5,17);
            %ia = ia(id);

            data = data(ia,:);
            cl = cell2mat(data(:,2));
            hi = cell2mat(data(:,4));
            lo = cell2mat(data(:,5));
            op = cell2mat(data(:,3));

            stdretC2C90d=backshift(1, smartMovingStd(calculateReturns(cl, 1), 90));

            %longs= op >= backshift(1, hi).*(1+entryZscore*stdretC2C90d);
            %shorts=op <= backshift(1, lo).*(1-entryZscore*stdretC2C90d);
            %改变了方向 数据相差10倍左右
            longs= op <= backshift(1, hi).*(1+entryZscore*stdretC2C90d);
            shorts=op >= backshift(1, lo).*(1-entryZscore*stdretC2C90d);


            positions=zeros(size(cl));

            positions(longs)=1;
            positions(shorts)=-1;

            ret=positions.*(op-cl)./op;
            ret(isnan(ret))=0;

            fprintf(1, '%s APR=%10.4f Sharpe=%4.2f\n', symbol, prod(1+ret).^(252/length(ret))-1, mean(ret)*sqrt(252)/std(ret));
            % APR=    0.1327 Sharpe=1.44
            cumret=cumprod(1+ret)-1; % compounded ROE

            y_re = cumret;
            %setfigure
            t_str = sprintf('%s',key_str);
            h = figure_S53(y_re,tref,t_str);
            re = {y_re,t_str};
        end
    end
    
end

function update_A71data()
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
end