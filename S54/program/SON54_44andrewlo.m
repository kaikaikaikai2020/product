%生成csi300/日本225/香港每日的交易明细文件，覆盖已有的来减少存储空间
clear;
close all
key_str = 'A44andrewlo综合计算';
path_save='计算结果';
if ~exist(path_save,'dir')
    mkdir(path_save)
end

addpath(genpath(fullfile(pwd,'jplv7')))
load A41Bog_data.mat

ind_sel = [1,2,4,6];
index=index(ind_sel);
index_0=index_0(ind_sel);
index_com=index_com(ind_sel);
index_db=index_db(ind_sel);
index_info = index_info(ind_sel);

index_pool = index;
index1 = cellfun(@(x) ['I_',x],index,'UniformOutput',false);
T = length(index_pool);
%T = 2;

re2 = cell(T,2);
re3 = re2;

for index_id = 1:T
    index_sel = index_pool{index_id};
    stocks = index_com{index_id};
    dtype=index_info{index_id};
    if strcmpi(dtype,'HK')
        fee1 = 1.2/1000;
        fee2 =1.2/1000;
    elseif strcmpi(dtype,'japan')
        fee1 = 1/10000;
        fee2 = 1/10000;
    elseif strcmpi(dtype,'us')
        fee1 = 1/10000;
        fee2 = 1/10000;
    elseif strcmpi(dtype,'csi')
        fee1 = 1/1000;
        fee2 = 1/1000;
    else
        sprintf('%s 未知品种',key_str)
        fee1 = 1/1000;
        fee2 = 1/1000;
    end    
    stks = get_pool_data(stocks,'2005-01-01',dtype);
    stocks=stks.stocks;
    %load tempstks.mat
    tday = stks.tday;
    cl=stks.cl;
    op=stks.op;
    tref = stks.tref;
    if size(tref,2)>size(tref,1)
        tref = tref';
    end
    ind = tday>=20100101;
    %ind = tday>=20070103 & tday<=20111230;
    tday=tday(ind);
    cl=cl(ind, :);
    op=op(ind, :);
    tref=tref(ind);
    % cl is a TxN array of closing prices, where T is the number of trading
    % days, and N is the number of stocks in the S&P 500
    ret=(cl-lag(cl, 1))./lag(cl, 1); % daily returns
    ret2 = (cl*(1-fee1)-lag(cl, 1)*(1+fee2))./(lag(cl, 1)*1+fee2); % daily returns real
    ret(ret>2)=0;
    ret2(ret2>0)=0;
    marketRet=smartmean(ret, 2); % equal weighted market index return

    weights=-(ret-repmat(marketRet, [1 size(ret, 2)]));
    %weights(weights<0)=0;
    weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

    %dailyret=smartsum(backshift(1, weights).*ret, 2); % Capital is always one
    dailyret=smartsum(backshift(1, weights).*ret2, 2); % Capital is always one
    dailyret(isnan(dailyret))=0;
    
    n=size(cl,2);
    info3 = cell(n,1);
    for i = 1:n
        temp = [tref,tref,tref,tref,num2cell([lag(cl(:,i), 1),cl(:,i),backshift(1, weights(:,i))])];
        temp(:,1)  = {'p1'};
        temp(:,2)={index_sel};
        temp(:,3) = stocks(i);
        ind= ~isnan(lag(cl(:,i),1)) & ~eq(lag(cl(:,i),1),0);
        temp = temp(ind,:);
        info3{i} = temp';
    end
    info3 = [info3{:}];
    
    %plot(cumprod(1+dailyret)-1); % Cumulative compounded return
    if ~isempty(tref)
        y_re = cumprod(1+dailyret);
        %setfigure
        figure
        subplot(2,1,1)
        h = figure_S53(y_re,tref,sprintf('%s-%s-part1',key_str,index_sel),0);
    else
        figure;
        plot(cumprod(1+dailyret)-1); % Cumulative compounded return
    end
    re2{index_id,1} = y_re;
    fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+dailyret).^(252/length(dailyret))-1, sqrt(252)*mean(dailyret)/std(dailyret));
    
    % % switch to use open prices
    % 
    ret=(op-backshift(1, cl))./backshift(1, cl); % daily returns
    ret(ret>2)=0;
    %ret2=(op*(1-fee1)-backshift(1, cl)*(1+fee2))./(backshift(1, cl)*(1+fee2); % daily returns

    marketRet=smartmean(ret, 2); % equal weighted market index return

    weights=-(ret-repmat(marketRet, [1 size(ret, 2)])); % weight of a stock is proportional to the negative distance to the market index.
    %weights(weights<0)=0;
    weights=weights./repmat(smartsum(abs(weights), 2), [1 size(weights, 2)]);

    %dailyret=smartsum(weights.*(cl-op)./op, 2)./smartsum(abs(weights), 2);
    dailyret=smartsum(weights.*(cl*(1-fee1)-op*(1+fee2))./(op*(1+fee2)), 2)./smartsum(abs(weights), 2);
    
    n=size(cl,2);
    info4 = cell(n,1);
    for i = 1:n
        temp = [tref,tref,tref,tref,num2cell([op(:,i),cl(:,i),weights(:,i)])];
        temp(:,1)  = {'p2'};
        temp(:,2)={index_sel};
        temp(:,3) = stocks(i);
        ind= ~isnan(op(:,i));
        temp = temp(ind,:);
        info4{i} = temp';
    end
    info4 = [info4{:}];
    
    
    dailyret(isnan(dailyret))=0;

    %plot(cumprod(1+dailyret)-1); % Cumulative compounded return
    y_re = cumprod(1+dailyret);
    %setfigure
    subplot(2,1,2);
    h = figure_S53(y_re,tref,sprintf('%s-%s-part2',key_str,index_sel),0);

    re2{index_id,2} = y_re;
    fprintf(1, 'APR=%f Sharpe=%f\n', prod(1+dailyret).^(252/length(dailyret))-1, sqrt(252)*mean(dailyret)/std(dailyret));
    
    sub_info = cellfun(@(x) ['A',x],stks.stocks,'UniformOutput',false);
    
    %pos = [tref,num2cell(weights)];
    %pos = cell2table(pos,'VariableNames',['date',sub_info']);
    %writetable(pos,sprintf('A_44_andrewlo_HK%s.csv',index_sel))
    re3{index_id}=[info3,info4];
    
end

sta_re = curve_static_batch(re2(:,2),index_info);

re3 = [re3{:}]';
[~,ia] = sort(re3(:,4));
re3 = re3(ia,:);
re3 = cell2table(re3,'VariableNames',{'type','code','code2','date','buy','sel','w'});
re3(isnan(re3.w),:) = [];

fn00=fullfile(path_save,sprintf('%s-历史交易明细.csv',key_str));
writetable(re3,fn00)

%writetable(pos,sprintf('A_44_andrewlo_csi%s.csv',index_sel))