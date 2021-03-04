%{
1）观测窗口是整个一天的数据
2）假设尾盘以close价格买入
3）统计第二天的收益来确定买卖
4）以第二天close价格卖出
5）不设其他stop loss规则
6）手续费还是双边千分之1.5
%}
%% SMTTradingModel
% SMT时域分形策略
% 参数设置--手续费设置
% 数据处理--剔除股指期货的前后15分钟数据
% 交易模拟--选取相关性前M个匹配
function fractalPostion =SMTTradingModelTool_update1(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate,P)
    %% 参数设置
    matchRecord=P.matchRecord;%匹配数据源：沪深300
    tradeRecord=P.tradeRecord;%交易数据源：股指期货主力合约
    tradeMin=P.tradeMin;%使用早盘120分钟K线数据进行分形匹配
    M=P.M;%找M个最为相似的交易日
    muchPara=P.muchPara;%多数上涨或下跌比例
    deanMethod=P.deanMethod;%1相关系数/2欧式距离/3兰氏距离/4曼哈顿距离
    testStart=P.testStart;
    cut_return = P.cut_return;
    %endDate='2018-6-1';%可设置最后一个交易日，目前数据最新到2018-5-31
    %% 数据处理
    % 输出匹配数据源和交易数据源
    [~,matchClose,matchDate]=dataSourse(matchRecord,tradeRecord,hsOpen,hsClose,hsDate,ifOpen,ifClose,ifDate);
    %% 交易模拟
    % 匹配数据窗口准备
    outMatchIndex=find(floor(matchDate)>=max(datenum(testStart),floor(ifDate(1))));
    outMatchClose=matchClose(outMatchIndex(1):end);
    outMatchDate=matchDate(outMatchIndex(1):end);
    
    inMatchDate=matchDate(1:outMatchIndex(1)-1);
    allDay=[1;find(floor(matchDate(1:end-1))~=floor(matchDate(2:end)))+1];
    outDay=[1;find(floor(outMatchDate(1:end-1))~=floor(outMatchDate(2:end)))+1];
    inDay=[1;find(floor(inMatchDate(1:end-1))~=floor(inMatchDate(2:end)))+1];
    %相关性匹配
    %corrRecord{1,1}:1 日期/2 相关系数/3 多空仓(1为多,-1做空)/ 4 尾盘最大逆向涨跌幅  5当日中午后涨跌
    corrRecord=corrSimilarityMatch2(allDay,outDay,inDay,matchClose,outMatchClose,matchDate,tradeMin,M,deanMethod);
    % 计算信号
    fractalPostion=dynamicTrailingStop2(corrRecord,outDay,outMatchDate,outMatchClose,tradeMin,muchPara,cut_return);
    fractalPostion = fractalPostion(:,1:2);
end

function corrRecord=corrSimilarityMatch2(allDay,outDay,inDay,matchClose,outMatchClose,matchDate,tradeMin,M,deanMethod)
%% 函数说明
% 匹配历史数据的相关性，筛选出相关系数超过一定值的日期，并记录
% corrRecord{1,1}:1 日期/2 相关系数/3 多空仓(1为多,-1做空)/ 4 尾盘最大逆向涨跌幅  5当日中午后涨跌
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
        elseif deanMethod==5          
            % 5dynamic time warping 举例
            dtw_re=dtw(x,y);
            corrAll(j,1)=dtw_re;
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
        % 判断下日尾盘情况
        if matchClose(allDay(j)+tradeMin,1)-matchClose(allDay(j+1)+tradeMin-1,1)>=0
            % 做空
            corrRecord{i,1}(ind_i,3)=-1;            
        elseif matchClose(allDay(j)+tradeMin,1)-matchClose(allDay(j+1)+tradeMin-1,1)<0
            % 做多
            corrRecord{i,1}(ind_i,3)=1;            
        end
        corrRecord{i,1}(ind_i,5) = matchClose(allDay(j+1)+tradeMin-1,1)/matchClose(allDay(j)+tradeMin,1)-1;
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

function fractalPostion=dynamicTrailingStop2(corrRecord,outDay,outMatchDate,outMatchClose,tradeMin,muchPara,cut_return)
%% 函数说明
% 动态追踪止损:判断当前亏损幅度是否超过历史匹配分形的最大亏损幅度均值，如超过则认为分形出现偏差
% 输出fractalPostion:1交易日期/2交易方向/3开仓下标/4历史最大亏损幅度均值/5当天最大亏损/6平仓下标
ind_j=1;
fractalPostion(ind_j,1)=0;
for i=1:length(outDay)
    %disp(['开始计算样本外第', num2str(i),'天的仓位情况——' datestr(floor(outMatchDate(outDay(i))),'yyyy-mm-dd')]);
    lenCorr=size(corrRecord{i,1},1); 
    longIndex=find(corrRecord{i,1}(:,3)==1);
    %增加了一个判断
    if length(longIndex)/lenCorr>muchPara
        % 历史匹配分形--做多
        fractalPostion(ind_j,1)=floor(outMatchDate(outDay(i)));%交易日期
        fractalPostion(ind_j,2)=1;%交易方向
        fractalPostion(ind_j,3)=outDay(i)+tradeMin-1;%开仓下标
        ind_j=ind_j+1;
    else
        % 历史匹配分形--做空
        fractalPostion(ind_j,1)=floor(outMatchDate(outDay(i)));%交易日期
        fractalPostion(ind_j,2)=-1;%交易方向
        fractalPostion(ind_j,3)=outDay(i)+tradeMin-1;%开仓下标
        ind_j=ind_j+1;
    end
end
end


