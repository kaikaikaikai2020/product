%% 封装原来程序，便于后续接入数据库
%修改下画图程序
%升级导出以下权重
%每次算出来的概率作为权重，比如计算出来0.6的可能性涨，我们就以0.6作为权重。。。
%如果下一个股票0.55的概率涨，权重就为0.55.
%升级原来框架，原来框架少算了一天
%% SMTTradingModel
% SMT时域分形策略
% 参数设置--手续费设置
% 数据处理--剔除股指期货的前后15分钟数据
% 交易模拟--选取相关性前M个匹配
function [tradeYield,result,tradeDetail,yearDetail,h,assure_ratio] =SMTTradingModelTool2(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate,P)
    %% 参数设置
    feeOpen=P.feeOpen;
    feeClose=P.feeClose;
    matchRecord=P.matchRecord;%匹配数据源：沪深300
    tradeRecord=P.tradeRecord;%交易数据源：股指期货主力合约
    tradeMin=P.tradeMin;%使用早盘120分钟K线数据进行分形匹配
    dayMin=P.dayMin;%每个交易日共240根1分钟K线
    M=P.M;%找M个最为相似的交易日
    muchPara=P.muchPara;%多数上涨或下跌比例
    deanMethod=P.deanMethod;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离
    stopMethod=P.stopMethod;%1收盘价止损/2触发价止损，不考虑能否交易/3跳开则开盘价止损，否则触发价止损
    testStart=P.testStart;
    trade_mode = P.trade_mode;
    if isfield(P,'out_sel')
        out_sel = P.out_sel;
    else
        out_sel = true;
    end
    %endDate='2018-6-1';%可设置最后一个交易日，目前数据最新到2018-5-31
    %% 数据处理
    % 输出匹配数据源和交易数据源
    [matchOpen,matchClose,matchDate,tradeOpen,tradeClose,tradeDate]=dataSourse(matchRecord,tradeRecord,hsOpen,hsClose,hsDate,ifOpen,ifClose,ifDate);
    %% 交易模拟
    % 匹配数据窗口准备
    outMatchIndex=find(floor(matchDate)>=max(datenum(testStart),floor(ifDate(1))));
    outMatchClose=matchClose(outMatchIndex(1):end);
    outMatchOpen=matchOpen(outMatchIndex(1):end);
    outMatchDate=matchDate(outMatchIndex(1):end);
    
    inMatchDate=matchDate(1:outMatchIndex(1)-1);
    allDay=[1;find(floor(matchDate(1:end-1))~=floor(matchDate(2:end)))+1];
    outDay=[1;find(floor(outMatchDate(1:end-1))~=floor(outMatchDate(2:end)))+1];
    inDay=[1;find(floor(inMatchDate(1:end-1))~=floor(inMatchDate(2:end)))+1];
    % 交易数据窗口准备
    outTradeIndex=find(floor(tradeDate)>=datenum(testStart));
    outTradeClose=tradeClose(outTradeIndex(1):end);
    outTradeOpen=tradeOpen(outTradeIndex(1):end);
    outTradeDate=tradeDate(outTradeIndex(1):end);

    % 相关性匹配
    corrRecord=corrSimilarityMatch2(allDay,outDay,inDay,matchClose,outMatchClose,matchDate,outMatchDate,tradeMin,M,deanMethod);
    % 动态止损交易
    [fractalPostion,assure_ratio]=dynamicTrailingStop2(corrRecord,outDay,outMatchDate,outMatchClose,tradeMin,muchPara);

    % 收益计算
    if stopMethod==1
        % 1、跌破止损价的第一根K线的收盘价
        tradeYield=tradeSimulate1(fractalPostion,outTradeClose,outTradeDate,feeOpen,feeClose);
    elseif stopMethod==2
        % 2、触发止损则平仓，不考虑价格是否能够交易
        tradeYield=tradeSimulate2(fractalPostion,outTradeClose,outTradeDate,feeOpen,feeClose,dayMin);
    elseif stopMethod==3
        % 3、如果平仓K线的开盘已跌破止损价，则使用开盘价止损，否则使用当根K线刚跌破止损价的点位止损
        tradeYield=tradeSimulate3(fractalPostion,outTradeClose,outTradeDate,feeOpen,feeClose,dayMin,outTradeOpen);
    end
    if out_sel
        % 结果输出
        [result,h]=outResultCaculate(tradeYield,feeOpen,feeClose,trade_mode);
        [tradeDetail,yearDetail]=outputSMT(tradeYield,outTradeDate);
    else
        result = [];
        h = [];
        tradeDetail = [];
        yearDetail = [];
    end
end

function corrRecord=corrSimilarityMatch2(allDay,outDay,inDay,matchClose,outMatchClose,matchDate,outMatchDate,tradeMin,M,deanMethod)
%% 函数说明
% 匹配历史数据的相关性，筛选出相关系数超过一定值的日期，并记录
% corrRecord{1,1}:日期/相关系数/多空仓(1为多,-1做空)/尾盘最大逆向涨跌幅
%%
corrRecord =cell(length(outDay),1);
for i=1:length(outDay)
    %disp(['开始匹配样本外第', num2str(i),'天数据——' datestr(floor(outMatchDate(outDay(i))),'yyyy-mm-dd')]);
    corrAll=zeros(length(inDay)+i-1,1);
    for j=1:length(inDay)+i-1
        x=outMatchClose(outDay(i):outDay(i)+tradeMin-1,1)/outMatchClose(outDay(i),1);
        y=matchClose(allDay(j):allDay(j)+tradeMin-1,1)/matchClose(allDay(j),1);
        if deanMethod==1
            % 1相关系数
            corr=corrcoef(x,y);
            corrAll(j,1)=corr(1,2);
        elseif deanMethod==2
            % 2欧氏距离
            euclidean=sum((x-y).^2);
            corrAll(j,1)=euclidean;
        elseif deanMethod==3            
            % 3 兰氏距离
            lsdean=sum(abs(x-y)./(x+y));
            corrAll(j,1)=lsdean;
        elseif deanMethod==4          
            % 4曼哈顿距离
            mhatondean=sum(abs(x-y));
            corrAll(j,1)=mhatondean;
        end
    end
    if deanMethod==1
        [~,index]=sort(corrAll,'descend');
    else
        [~,index]=sort(corrAll,'ascend');
    end
    for ind_i=1:M
        j=index(ind_i);
        corrRecord{i,1}(ind_i,1)=floor(matchDate(allDay(j)));
        corrRecord{i,1}(ind_i,2)=corrAll(j,1);
        % 判断当日尾盘情况
        if matchClose(allDay(j)+tradeMin-1,1)-matchClose(allDay(j+1)-1,1)>=0
            % 做空则记录最大涨幅
            corrRecord{i,1}(ind_i,3)=-1;
            corrRecord{i,1}(ind_i,4)=max(matchClose(allDay(j)+tradeMin-1:allDay(j+1)-1,1))/matchClose(allDay(j)+tradeMin-1,1)-1;
        elseif matchClose(allDay(j)+tradeMin-1,1)-matchClose(allDay(j+1)-1,1)<0
            % 做多则记录最大跌幅
            corrRecord{i,1}(ind_i,3)=1;
            corrRecord{i,1}(ind_i,4)=1-min(matchClose(allDay(j)+tradeMin-1:allDay(j+1)-1,1))/matchClose(allDay(j)+tradeMin-1,1);
        end
    end
end
end
function [matchOpen,matchClose,matchDate,tradeOpen,tradeClose,tradeDate]=dataSourse(matchRecord,tradeRecord,hsOpen,hsClose,hsDate,ifOpen,ifClose,ifDate)
%% 函数说明
% 根据输入参数，判断沪深300和股指期货的匹配数据源和交易数据源
% 两部分数据，一部分用于计算匹配度，一部分用于回测（一部分发信号，一部分执行交易，两部分可以相同，也可以不同）
%%
% 判断匹配数据源
if matchRecord==1
    matchOpen=hsOpen;
    matchClose=hsClose;
    matchDate=hsDate;
else
    matchOpen=ifOpen;
    matchClose=ifClose;
    matchDate=ifDate;
end
% 判断交易数据源
if tradeRecord==1
    tradeOpen=hsOpen;
    tradeClose=hsClose;
    tradeDate=hsDate;
else
    tradeOpen=ifOpen;
    tradeClose=ifClose;
    tradeDate=ifDate;
end
end
function dateVec=dateStandard(dateSeries)
%% 函数说明
% 将时间序列标准化
% 输入数值型的时间序列，输出为时间向量
%% 
dateVec=datevec(dateSeries);
dateVec(:,5)=dateVec(:,5)+round(dateVec(:,6)/60);
dateVec(:,6)=0;

end

function [fractalPostion,assure_ratio]=dynamicTrailingStop2(corrRecord,outDay,outMatchDate,outMatchClose,tradeMin,muchPara)
%% 函数说明
% 动态追踪止损:判断当前亏损幅度是否超过历史匹配分形的最大亏损幅度均值，如超过则认为分形出现偏差
% 输出fractalPostion:1交易日期/2交易方向/3开仓下标/4历史最大亏损幅度均值/5当天最大亏损/6平仓下标
%assure_ratio 把握度，上涨或下跌的概率
ind_j=1;
fractalPostion(ind_j,1)=0;
assure_ratio = ones(length(outDay)-1,1);
for i=1:length(outDay)
    %disp(['开始计算样本外第', num2str(i),'天的仓位情况——' datestr(floor(outMatchDate(outDay(i))),'yyyy-mm-dd')]);
    lenCorr=size(corrRecord{i,1},1); 
    longIndex=find(corrRecord{i,1}(:,3)==1);
    shortIndex=find(corrRecord{i,1}(:,3)==-1);
    if length(longIndex)/lenCorr>muchPara
        % 历史匹配分形--做多
        fractalPostion(ind_j,1)=floor(outMatchDate(outDay(i)));%交易日期
        fractalPostion(ind_j,2)=1;%交易方向
        fractalPostion(ind_j,3)=outDay(i)+tradeMin-1;%开仓下标
        lossMeanRate=mean(corrRecord{i,1}(longIndex,4));            
        fractalPostion(ind_j,4)=lossMeanRate;%历史最大扣亏损幅度均值
        if i < length(outDay)
            lossNowRate=1-outMatchClose(outDay(i)+tradeMin:outDay(i+1)-1)./outMatchClose(outDay(i)+tradeMin-1);
        else
            lossNowRate=1-outMatchClose(outDay(i)+tradeMin:end)./outMatchClose(outDay(i)+tradeMin-1);
        end
        indC=find(lossNowRate>=lossMeanRate);
        fractalPostion(ind_j,5)=max(lossNowRate);%当天最大亏损幅度
        if isempty(indC)==1
            if i < length(outDay)
                fractalPostion(ind_j,6)=outDay(i+1)-1;%平仓下标
            else
                fractalPostion(ind_j,6)=length(outMatchClose);%平仓下标
            end
        else
            fractalPostion(ind_j,6)=indC(1)+outDay(i)+tradeMin-1;%平仓下标
        end
        assure_ratio(ind_j,1) = length(longIndex)/(length(longIndex)+length(shortIndex));
        ind_j=ind_j+1;
    else
        % 历史匹配分形--做空
        fractalPostion(ind_j,1)=floor(outMatchDate(outDay(i)));%交易日期
        fractalPostion(ind_j,2)=-1;%交易方向
        fractalPostion(ind_j,3)=outDay(i)+tradeMin-1;%开仓下标
        lossMeanRate=mean(corrRecord{i,1}(shortIndex,4));            
        fractalPostion(ind_j,4)=lossMeanRate;%历史最大扣亏损幅度均值
        if i < length(outDay)
            lossNowRate=outMatchClose(outDay(i)+tradeMin:outDay(i+1)-1)./outMatchClose(outDay(i)+tradeMin-1)-1;
        else
            lossNowRate=outMatchClose(outDay(i)+tradeMin:end)./outMatchClose(outDay(i)+tradeMin-1)-1;
        end
        indC=find(lossNowRate>=lossMeanRate);
        fractalPostion(ind_j,5)=max(lossNowRate);%当天最大亏损幅度
        if isempty(indC)==1
            if i < length(outDay)
                fractalPostion(ind_j,6)=outDay(i+1)-1;%平仓下标
            else
                fractalPostion(ind_j,6)=length(outMatchClose);%平仓下标
            end
        else
            fractalPostion(ind_j,6)=indC(1)+outDay(i)+tradeMin-1;%平仓下标
        end
        assure_ratio(ind_j,1) = length(shortIndex)/(length(longIndex)+length(shortIndex));
        ind_j=ind_j+1;
    end
end
end

function [tradeDetail,yearDetail]=outputSMT(tradeYield,outTradeDate)
%% 函数说明
% 输出每笔交易详情和历史年度收益到Excel
len=size(tradeYield,1);
tradeDetail{len,1}=[];
for i=1:len
    tradeDetail{i,1}=datestr(tradeYield(i,1),'yyyy-mm-dd');
    tradeDetail{i,2}=datestr(outTradeDate(tradeYield(i,6)),'HH:MM');
    tradeDetail{i,5}=datestr(outTradeDate(tradeYield(i,7)),'HH:MM');
    tradeDetail{i,3}=tradeYield(i,4);
    tradeDetail{i,6}=tradeYield(i,5);    
    if tradeYield(i,3)==1
        tradeDetail{i,4}='多';
    else
        tradeDetail{i,4}='空';
    end
    tradeDetail{i,7}=[num2str(floor(tradeYield(i,2)*10000)/100),'%'];
end
tradeDetail=[{'交易日期','开仓时间','开仓点位','开仓方向','平仓时间','平仓点位','收益率'};tradeDetail];
hisYear=unique(year(tradeYield(:,1)));
yearDetail{1,1}=[];
for i=1:length(hisYear)
    index=find(year(tradeYield(:,1))==hisYear(i));
    yield=cumprod(tradeYield(index,2)+1)-1;
    yearDetail{i,1}=hisYear(i);
    yearDetail{i,2}=length(find(tradeYield(index,3)~=0));
    yearDetail{i,3}=[num2str(floor(yield(end)*10000)/100),'%'];
    yearDetail{i,4}=[num2str(floor(sum(tradeYield(index,2)>0)/yearDetail{i,2}*10000)/100),'%'];    
    yearDetail{i,5}=[num2str(floor(maxdrawdown(yield+1)*10000)/100),'%'];
    yearDetail{i,6}=sum(tradeYield(index,2)>0);
    yearDetail{i,7}=sum(tradeYield(index,2)<0);
end
yearDetail=[{'年度','交易次数','累计收益率','成功率','最大回撤','盈利次数','亏损次数'};yearDetail];
end
function [result1,h]=outResultCaculate(tradeYield,feeOpen,feeClose,trade_mode)
%% 函数说明
% 输出当前交易结果及统计数据
%%
allDayNum=length(tradeYield(:,1));
longNum=length(find(tradeYield(:,3)==1));
shortNum=length(find(tradeYield(:,3)==-1));
tradeNum=longNum+shortNum;
winNum=length(find(tradeYield(:,2)>0));
loseNum=length(find(tradeYield(:,2)<0));
winRate=mean(tradeYield(tradeYield(:,2)>0,2));
loseRate=mean(tradeYield(tradeYield(:,2)<0,2));
%只多仓 修正
if eq(trade_mode,1)
    Yield=tradeYield(:,2);
    Yield(tradeYield(:,3)<0)=0;
    % 不考虑交易费用的累计收益率
    noYield=Yield;
    noYield(tradeYield(:,3)>0)=tradeYield(tradeYield(:,3)>0,2)+feeOpen+feeClose;
else
    Yield=tradeYield(:,2);
    noYield=Yield;
    noYield(abs(tradeYield(:,3))>0)=tradeYield(abs(tradeYield(:,3))>0,2)+feeOpen+feeClose;
end

cumYield=cumprod(Yield+1)-1;
yearYield=cumYield(end)/allDayNum*250;
maxDrawDown=maxdrawdown(cumYield+1);

noCumYield=cumprod(noYield+1)-1;
result=[allDayNum,tradeNum,winNum,loseNum,winRate,loseRate,winNum/tradeNum,cumYield(end),maxDrawDown,yearYield,noCumYield(end)];

result_info = {'总计交易日','交易次数','盈利次数','亏损次数','平均盈利率%','平均亏损率%',...
    '成功率%','累计收益率%','最大回撤率%','年化收益率%','累计收益（无手续费）%'};
result(5:end) = result(5:end)*100;
result1 = [result_info;num2cell(result)]';

%% 结果输出
% disp('结果输出如下：')
% disp(['总计交易日：' num2str(allDayNum)])
% disp(['交易次数：' num2str(tradeNum)])
% disp(['盈利次数：' num2str(winNum)])
% disp(['亏损次数：' num2str(loseNum)])
% disp(['平均盈利率：' num2str(floor(winRate*10000)/100) '%'])
% disp(['平均亏损率：' num2str(floor(loseRate*10000)/100) '%'])
% disp(['成功率：' num2str(floor(winNum/tradeNum*10000)/100) '%'])
% disp(['累计收益率：' num2str(floor(cumYield(end)*10000)/100) '%'])
% disp(['最大回撤率：' num2str(floor(maxDrawDown*10000)/100) '%'])
% disp(['年化收益率：' num2str(floor(yearYield*10000)/100) '%'])

%% 画图
%h = [];
h  = figure;
subplot(2,1,1)
plot(cumYield,'r')
hold on
plot(noCumYield,'b')
t_str = cellstr(datestr(tradeYield(:,1),'yyyymmdd'));
T = length(t_str);
set(gca,'xlim',[0,T]);
set(gca,'XTick',floor(linspace(1,T,10)));
set(gca,'XTickLabel',t_str(floor(linspace(1,T,10))));
set(gca,'XTickLabelRotation',90)
legend('有交易费','无交易费','location','best')
title({['累计收益' num2str(floor(cumYield(end)*10000)/100) '%，年化收益' num2str(floor(yearYield*10000)/100) '%，' '最大回撤' num2str(floor(maxDrawDown*10000)/100) '%' ],...
    ['总交易次数' num2str(tradeNum) '，平均盈利率' num2str(floor(winRate*10000)/100) '%，' '平均亏损率' num2str(floor(loseRate*10000)/100) '%']})

subplot(2,1,2)
%Yield
yearS_all = year(tradeYield(:,1));
yearS = unique(yearS_all);
yearRate = zeros(size(yearS));
for i = 1:length(yearS)
    sub_ind = eq(yearS_all,yearS(i));
    temp = cumprod(1+Yield(sub_ind))-1;
    yearRate(i) = temp(end);
end

redRate=yearRate;
redRate(redRate<=0)=0;
blueRate=yearRate;
blueRate(blueRate>0)=0;
bar(yearS,redRate,'r')
hold on
bar(yearS,blueRate,'b')
set(gca,'XTickLabelRotation',90)
title('有交易费用年度收益情况')

end
function tradeYield=tradeSimulate1(fractalPostion,outTradeClose,outTradeDate,feeOpen,feeClose)
%% 函数说明
% 根据给出的仓位对股指进行交易
% 止损操作——使用跌破止损价的第一根K线的收盘价进行止损操作
% 输出结果tradeYield:日期/当日收益/当日开仓方向/开仓点位/平仓点位/开仓时间/平仓时间
%%
outTradeDay=unique(floor(outTradeDate));
tradeNum=length(outTradeDay);
tradeYield(1,1)=0;
ind=find(outTradeDay==fractalPostion(1,1));
for i=ind:tradeNum
    index=find(fractalPostion(:,1)==outTradeDay(i,1));
    if isempty(index)==0
        tradeYield(i,1)=fractalPostion(index(1),2)*(outTradeClose(fractalPostion(index(1),6))/outTradeClose(fractalPostion(index(1),3))-1)-feeOpen-feeClose;
        tradeYield(i,2)=fractalPostion(index(1),2);
        tradeYield(i,3)=outTradeClose(fractalPostion(index(1),3));
        tradeYield(i,4)=outTradeClose(fractalPostion(index(1),6));
        tradeYield(i,5)=fractalPostion(index(1),3);
        tradeYield(i,6)=fractalPostion(index(1),6);
    else
        tradeYield(i,1)=0;
        tradeYield(i,2)=0;
    end
end
tradeYield=[outTradeDay,tradeYield];
tradeYield=tradeYield(ind:end-1,:);
end
function tradeYield=tradeSimulate2(fractalPostion,outTradeClose,outTradeDate,feeOpen,feeClose,dayMin)
%% 函数说明
% 根据给出的仓位对股指进行交易
% 止损操作——触发止损则平仓，不考虑价格是否能够交易
% 输出结果tradeYield:日期/当日收益/当日开仓方向/开仓点位/平仓点位/开仓时间/平仓时间
%%
outTradeDay=unique(floor(outTradeDate));
tradeNum=length(outTradeDay);
tradeYield(1,1)=0;
ind=find(outTradeDay==fractalPostion(1,1));
for i=ind:tradeNum 
    index=find(fractalPostion(:,1)==outTradeDay(i,1));
    if isempty(index)==0
        if mod(fractalPostion(index(1),6),dayMin)==0
            tradeYield(i,1)=fractalPostion(index(1),2)*(outTradeClose(fractalPostion(index(1),6))/outTradeClose(fractalPostion(index(1),3))-1)-feeOpen-feeClose;
            tradeYield(i,2)=fractalPostion(index(1),2);
            tradeYield(i,3)=outTradeClose(fractalPostion(index(1),3));
            tradeYield(i,4)=outTradeClose(fractalPostion(index(1),6));
            tradeYield(i,5)=fractalPostion(index(1),3);
            tradeYield(i,6)=fractalPostion(index(1),6);
        else            
            tradeYield(i,1)=-fractalPostion(index(1),4)-feeOpen-feeClose;
            tradeYield(i,2)=fractalPostion(index(1),2);
            tradeYield(i,3)=outTradeClose(fractalPostion(index(1),3));
            tradeYield(i,4)=tradeYield(i,3)*(1-tradeYield(i,2)*fractalPostion(index(1),4));    
            tradeYield(i,5)=fractalPostion(index(1),3);
            tradeYield(i,6)=fractalPostion(index(1),6);
        end
    else
        tradeYield(i,1)=0;
        tradeYield(i,2)=0;
    end
end
tradeYield=[outTradeDay,tradeYield];
tradeYield=tradeYield(ind:end-1,:);
end
function tradeYield=tradeSimulate3(fractalPostion,outTradeClose,outTradeDate,feeOpen,feeClose,dayMin,outTradeOpen)
%% 函数说明
% 根据给出的仓位对股指进行交易
% 止损操作——如果平仓K线的开盘已跌破止损价，则使用开盘价止损，否则使用当根K线刚跌破止损价的点位止损
% 输出结果tradeYield:日期/当日收益/当日开仓方向/开仓点位/平仓点位/开仓时间/平仓时间
%%
outTradeDay=unique(floor(outTradeDate));
tradeNum=length(outTradeDay);
tradeYield(1,1)=0;
ind=find(outTradeDay==fractalPostion(1,1));
%手续费统计方法升级
for i=ind:tradeNum
    index=find(fractalPostion(:,1)==outTradeDay(i,1));
    if isempty(index)==0
        if mod(fractalPostion(index(1),6),dayMin)==0
            tradeYield(i,1)=fractalPostion(index(1),2)*(outTradeClose(fractalPostion(index(1),6))/outTradeClose(fractalPostion(index(1),3))-1)-feeOpen-feeClose;
            tradeYield(i,2)=fractalPostion(index(1),2);
            tradeYield(i,3)=outTradeClose(fractalPostion(index(1),3));
            tradeYield(i,4)=outTradeClose(fractalPostion(index(1),6));
            tradeYield(i,5)=fractalPostion(index(1),3);
            tradeYield(i,6)=fractalPostion(index(1),6);
        else
            openRate=fractalPostion(index(1),2)*(outTradeOpen(fractalPostion(index(1),6))/outTradeClose(fractalPostion(index(1),3))-1);
            if -openRate>fractalPostion(index(1),4)
                tradeYield(i,1)=openRate-feeOpen-feeClose;
                tradeYield(i,2)=fractalPostion(index(1),2);
                tradeYield(i,3)=outTradeClose(fractalPostion(index(1),3));
                tradeYield(i,4)=outTradeOpen(fractalPostion(index(1),6));
                tradeYield(i,5)=fractalPostion(index(1),3);
                tradeYield(i,6)=fractalPostion(index(1),6);
            else
                tradeYield(i,1)=-fractalPostion(index(1),4)-feeOpen-feeClose;
                tradeYield(i,2)=fractalPostion(index(1),2);
                tradeYield(i,3)=outTradeClose(fractalPostion(index(1),3));
                tradeYield(i,4)=tradeYield(i,3)*(1-tradeYield(i,2)*fractalPostion(index(1),4));
                tradeYield(i,5)=fractalPostion(index(1),3);
                tradeYield(i,6)=fractalPostion(index(1),6);
            end
        end
    else
        tradeYield(i,1)=0;
        tradeYield(i,2)=0;
    end
end
tradeYield=[outTradeDay,tradeYield];
%tradeYield=tradeYield(ind:end-1,:);
end
