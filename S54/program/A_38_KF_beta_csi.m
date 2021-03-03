clear;
% Daily data on EWA-EWC
addpath(genpath(fullfile(pwd,'jplv7')))
key_str = 'A-38';
index_pool = {'000016','000300';'000300','000905';'000016','000905'};
index_info1={'I50','I300';'I300','I500';'I50','I500'};
index_info = {'50-300','300-500','50-500'};

T = size(index_pool,1);
re = cell(T,1);
tref0 = [];
for index_ID = 1:3
    sql_str = 'select tradeDate,closeIndex from yuqerdata.yq_index where symbol = "%s" order by tradeDate';
    
    sub_ticker1 = index_pool{index_ID,1};
    sub_ticker2 = index_pool{index_ID,2};
    sym1 = index_info1{index_ID,1};
    sym2 = index_info1{index_ID,2};
    
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
    
    %{
    longsEntry=e < -sqrt(Q); % a long position means we should buy EWC
    longsExit=e > -sqrt(Q);

    shortsEntry=e > sqrt(Q);
    shortsExit=e < sqrt(Q);
    %}
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

    if ~isempty(tref)
        y_re = cumprod(1+ret)-1;
        %setfigure
        h = figure_S53(y_re,tref,[]);
        %title(sprintf('%s-%s-%s',key_str,sub_ticker1,sub_ticker2))
        title(sprintf('%s-%s',key_str,index_info{index_ID}))
    else
        figure;
        plot(cumprod(1+ret)-1); % Cumulative compounded return
    end

    fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+ret).^(252/length(ret))-1, sqrt(252)*mean(ret)/std(ret));
    % APR=0.262252 Sharpe=2.361162

    %title('A-38-usdtwd vs usdkrw')
    
    
    pos0=[tref,num2cell([numUnits,positions])];
    pos0=cell2table(pos0,'VariableNames',{'date','signal',sym1,sym2});
    re{index_ID}=pos0;
    tref0=unique([tref0;table2cell(pos0(:,1))]);
    
end

X = cell(T,1);
var =cell(T,1);
for i = 1:T
    sub_x = re{i};
    temp = sub_x.Properties.VariableNames(2:end);
    temp = cellfun(@(x) sprintf('%s_I%d',x,i),temp,'UniformOutput',false);
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

writetable(X,'A_38_KF_beta_csi.csv')
