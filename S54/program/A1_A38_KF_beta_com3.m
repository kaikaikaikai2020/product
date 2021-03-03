clear
%所有结果整合
%综合程序
addpath(genpath(fullfile(pwd,'jplv7')))
%updateA38data()
load bloombergdata

tasklist = [1,2,3];
re = cell(3,1);
if any(eq(tasklist,1))
    figure
    re{1} = update_part1(data);
end
if any(eq(tasklist,2))
    figure
    re{2} = update_part2(data);
end
if any(eq(tasklist,3))
    figure
    re{3}=update_part3(data);    
end

ind = cellfun(@isempty,re);
re(ind) = [];
for i = 1:length(re)
    if eq(i,1)
        tref = re{i}(:,1);
    else
        tref = unique([tref;re{i}(:,1)]);
    end
end

X = tref;
for i = 1:length(re)
    X=outerjoin(X,re{i},'key','date','MergeKeys' ,1);
    %X.Properties.VariableNames{1}='date';
end

writetable(X,sprintf('A1_A38_KF_综合_%s.csv',datestr(now,'yymmdd')))
%%%%%%%
function X = update_part3(data)
symbol_pool = {'HSI','HSCEI';'NKY','TPX';...
    'TAMSCI','TWSE'};
%delta_pool = [1e-5,1e-10,1e-4];
T = size(symbol_pool,1);
delta_pool = [1e-5,1e-10,1e-4];
dir = [-1,1,1];
re = cell(T,1);
tref0 = [];
for pool_id = 1:T
    sym1 = symbol_pool{pool_id,1};
    sym2 = symbol_pool{pool_id,2};

    x=data.(sym1);
    y=data.(sym2);

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
    y2=[x(:, 1) y];
    
    if eq(dir(pool_id),1)
        longsEntry=e < -sqrt(Q); % a long position means we should buy EWC
        longsExit=e > -sqrt(Q);

        shortsEntry=e > sqrt(Q);
        shortsExit=e < sqrt(Q);
    else
        longsEntry=e > -sqrt(Q); % a long position means we should buy EWC
        longsExit=e < -sqrt(Q);

        shortsEntry=e < sqrt(Q);
        shortsExit=e > sqrt(Q);
    end

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

    if ~isempty(tref)
        y_re = cumprod(1+ret)-1;
        %setfigure
        subplot(3,1,pool_id)
        figure_S53_update(y_re,tref,[]);
    else
        figure;
        plot(cumprod(1+ret)-1); % Cumulative compounded return
    end

    fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+ret).^(252/length(ret))-1, sqrt(252)*mean(ret)/std(ret));
    % APR=0.262252 Sharpe=2.361162
    sym1=replace(sym1,'Curncy','');
    sym2=replace(sym2,'Curncy','');
    temp = sprintf('A-38-%s vs %s',sym1,sym2);
    temp = strrep(temp,'_','-');
    if eq(dir(pool_id),-1)
        temp = sprintf('%s-sig-reverse',temp);
    end
    title(temp)
    sprintf('position of %s-%s',sym1,sym2)
    %pos0=[tref,num2cell([numUnits,positions])];
    %pos0=cell2table(pos0,'VariableNames',{'date',sprintf('signal',sym1,sym2),sym1,sym2});
    temp1 = [numUnits,positions,y2];
    temp1 = temp1([1:end,end],:);
    temp1 = lag(temp1,1);
    temp2 = y2([1:end,end],:);
    temp2(end,:) = nan;
    tref = tref([1:end,end],:);
    tref{end} = '下一个交易日';

    pos0=[tref,num2cell([temp1,temp2])];
    pos0=cell2table(pos0,'VariableNames',{'date','signal',['pos_',sym1],['pos_',sym2],['buy_',sym1],['buy_',sym2],['sell_',sym1],['sell_',sym2]});
    %writetable(pos0,sprintf('%s%s-%s.csv',sym1,sym2,tref{end}))
    re{pool_id}=pos0;
    tref0=unique([tref0;table2cell(pos0(:,1))]);
end

X = cell(T,1);
var =cell(T,1);
for i = 1:T
    sub_x = re{i};
    temp = sub_x.Properties.VariableNames(2:end);
    temp = cellfun(@(x) sprintf('%s_P3I%d',x,i),temp,'UniformOutput',false);
    var{i} = temp;
    sub_x = table2cell(sub_x);
    sub_x2 = nan(length(tref0),size(sub_x,2)-1);
    [~,ia,ib] = intersect(tref0,sub_x(:,1),'stable');
    sub_x2(ia,:) = cell2mat(sub_x(ib,2:end));
    X{i} = sub_x2;
    
end

X=[tref0,num2cell([X{:}])];
var=['date',[var{:}]];
X = cell2table(X,'VariableNames',var);
%writetable(X,sprintf('A1_A38_KF_beta_part3_%s.csv',datestr(now,'yymmddHHMM')))
%writetable(X,'A1_A38_KF_beta_part3.csv')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = update_part2(data)
symbol_pool = {'KWN1M','NTN1M';'NTN1M','NTN2M';...
    'NTN1M','NTN3M';'NTN1M','NTN6M';...
    'NTN1M','NTN12M';'NTN6M','NTN12M'};
%delta_pool = [1e-5,1e-10,1e-4];
T = size(symbol_pool,1);
delta_pool=[1e-8,1e-8,1e-8,1e-8,1e-8,1e-8];
re = cell(T,1);
tref0 = [];
for pool_id = 1:T
    sym1 = symbol_pool{pool_id,1};
    sym2 = symbol_pool{pool_id,2};

    x=data.(sym1);
    y=data.(sym2);

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
    
    temp = sprintf('A-38-%s vs %s',sym1,sym2);
    temp = strrep(temp,'_','-');
    title(temp)
    
    if ~isempty(tref)
        y_re = cumprod(1+ret)-1;
        %setfigure
        subplot(3,2,pool_id)
        figure_S53_update(y_re,tref,temp);
    else
        figure;
        plot(cumprod(1+ret)-1); % Cumulative compounded return
    end

    fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+ret).^(252/length(ret))-1, sqrt(252)*mean(ret)/std(ret));
    % APR=0.262252 Sharpe=2.361162
    sym1=replace(sym1,'Curncy','');
    sym2=replace(sym2,'Curncy','');
    
    sprintf('position of %s-%s',sym1,sym2)
    %pos0=[tref,num2cell([numUnits,positions])];
    %pos0=cell2table(pos0,'VariableNames',{'date',sprintf('signal',sym1,sym2),sym1,sym2});
    temp1 = [numUnits,positions,y2];
    temp1 = temp1([1:end,end],:);
    temp1 = lag(temp1,1);
    temp2 = y2([1:end,end],:);
    temp2(end,:) = nan;
    tref = tref([1:end,end],:);
    tref{end} = '下一个交易日';

    pos0=[tref,num2cell([temp1,temp2])];
    pos0=cell2table(pos0,'VariableNames',{'date','signal',['pos_',sym1],['pos_',sym2],['buy_',sym1],['buy_',sym2],['sell_',sym1],['sell_',sym2]});
    
    %writetable(pos0,sprintf('%s%s-%s.csv',sym1,sym2,tref{end}))
    re{pool_id}=pos0;
    tref0=unique([tref0;table2cell(pos0(:,1))]);
end

X = cell(T,1);
var =cell(T,1);
for i = 1:T
    sub_x = re{i};
    temp = sub_x.Properties.VariableNames(2:end);
    temp = cellfun(@(x) sprintf('%s_P2I%d',x,i),temp,'UniformOutput',false);
    var{i} = temp;
    sub_x = table2cell(sub_x);
    sub_x2 = nan(length(tref0),size(sub_x,2)-1);
    [~,ia,ib] = intersect(tref0,sub_x(:,1),'stable');
    sub_x2(ia,:) = cell2mat(sub_x(ib,2:end));
    X{i} = sub_x2;
    
end

X=[tref0,num2cell([X{:}])];
var=['date',[var{:}]];
X = cell2table(X,'VariableNames',var);
%writetable(X,sprintf('A1_A38_KF_beta_part2_%s.csv',datestr(now,'yymmddHHMM')))
%writetable(X,'A1_A38_KF_beta_part2.csv')
end


function X = update_part1(data)
    symbol_pool = {'NTN1M','CNH1M';'CNH1M','CNH2M';...
        'CNH1M','CNH3M';'CNH1M','CNH6M';...
        'CNH1M','CNH12M';'CNH6M','CNH12M';...
        'NTN1M','IRN1M';'IRN1M','IRN2M';...
        'IRN1M','IRN3M';'IRN1M','IRN6M';...
        'IRN1M','IRN12M';'IRN6M','IRN12M'};
    %delta_pool = [1e-5,1e-10,1e-4];
    T = size(symbol_pool,1);
    %delta_pool=[1e-8,1e-8,1e-8,1e-8,1e-8,1e-8];
    delta_pool= [1e-7,1e-10,1e-10,...
        1e-10,1e-10,1e-10,...
        ones(1,6).*1e-5];
    re = cell(T,1);
    tref0 = [];
    for pool_id = 1:T
        sym1 = symbol_pool{pool_id,1};
        sym2 = symbol_pool{pool_id,2};

        x=data.(sym1);
        y=data.(sym2);

        x(cellfun(@isnan,x(:,2)),:) = [];
        y(cellfun(@isnan,y(:,2)),:) = [];


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

        if ~isempty(tref)
            y_re = cumprod(1+ret)-1;
            %setfigure
            subplot(3,4,pool_id)
            %figure
            figure_S53_update(y_re,tref,[]);
        else
            figure;
            plot(cumprod(1+ret)-1); % Cumulative compounded return
        end
        sym1=replace(sym1,'Curncy','');
        sym2=replace(sym2,'Curncy','');

        fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+ret).^(252/length(ret))-1, sqrt(252)*mean(ret)/std(ret));
        % APR=0.262252 Sharpe=2.361162
        temp = sprintf('A-38-%s vs %s',sym1,sym2);
        temp = strrep(temp,'_','-');
        title(temp)
        sprintf('position of %s-%s',sym1,sym2)
        %pos0=[tref,num2cell([numUnits,positions,lag(y2),y2])];
        %pos0=pos0([1:end,end],:);
        temp1 = [numUnits,positions,y2];
        temp1 = temp1([1:end,end],:);
        temp1 = lag(temp1,1);
        temp2 = y2([1:end,end],:);
        temp2(end,:) = nan;
        tref = tref([1:end,end],:);
        tref{end} = '下一个交易日';
        
        pos0=[tref,num2cell([temp1,temp2])];
        pos0=cell2table(pos0,'VariableNames',{'date','signal',['pos_',sym1],['pos_',sym2],['buy_',sym1],['buy_',sym2],['sell_',sym1],['sell_',sym2]});
        %writetable(pos0,sprintf('%s%s-%s.csv',sym1,sym2,tref{end}))
        re{pool_id}=pos0;
        tref0=unique([tref0;table2cell(pos0(:,1))]);
    end

    X = cell(T,1);
    var =cell(T,1);
    for i = 1:T
        sub_x = re{i};
        temp = sub_x.Properties.VariableNames(2:end);
        temp = cellfun(@(x) sprintf('%s_P1I%d',x,i),temp,'UniformOutput',false);
        var{i} = temp;
        sub_x = table2cell(sub_x);
        sub_x2 = nan(length(tref0),size(sub_x,2)-1);
        [~,ia,ib] = intersect(tref0,sub_x(:,1),'stable');
        sub_x2(ia,:) = cell2mat(sub_x(ib,2:end));
        X{i} = sub_x2;    
    end

    X=[tref0,num2cell([X{:}])];
    var=['date',[var{:}]];
    X = cell2table(X,'VariableNames',var);
    %writetable(X,sprintf('A1_A38_KF_beta_part1_%s.csv',datestr(now,'yymmddHHMM')))
    %writetable(X,'A1_A38_KF_beta_part1.csv')
end