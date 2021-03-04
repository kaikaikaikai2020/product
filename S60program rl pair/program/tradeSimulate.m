function [accountDetail,tradeDetail]=tradeSimulate(position,priceH,priceA,dateS,V0,feeH,feeA,sRateH,sRateA,stampH,stampA,r,cDay)
%% 函数说明
% 根据计算的最优配置仓位，每日对股票对进行模拟交易，尽量模拟真实交易
% 初始投资金额为V0,融券利率算头不算尾按照自然天数计算
% 无风险利息算头不算尾，按照自然天数来计算
% accountDetail:1总资产/2H市值/3H持仓数量及方向/4H交易费用/5H融券利息/6A市值/7A持仓数量及方向/8A交易费用/9A融券利息/10可用资金/11无风险利息/12H融券额度/13A融券额度
% tradeDetail:1H交易数量及方向/2H交易价格/3A交易数量及方向/4A交易价格/5H交易盈利/6A交易盈利
%%
perH=100;
perA=100;
shortMonth=6;%预留6个月融券费用（融券最长期限为6个月）
totalT=length(dateS);
startT=min(find(position~=0));
totalAsset(1:startT-1,1)=V0;
accountDetail(1:startT-1,13)=0;
accountDetail(1:startT-1,[1,10])=V0;%初期未建仓前不算无风险利息
tradeDetail(1:startT-1,1:6)=0;

% 初始交易
if position(startT)>0
    %买开H股，融券卖开A股
    tradeFlagH=1;
    tradeFlagA=-1;
    preValue=position(startT)*totalAsset(startT-1,1)/(1+feeH+feeA+stampA+sRateA*shortMonth/12);%预留手续费(买开H，卖开A)+shortMonth个月的融券费用
    sharesH=floor(preValue/priceH(startT,1)/perH)*perH;
    marketH=sharesH*priceH(startT,1);
    tradeFeeH=marketH*feeH;
    shortFeeH=0;
    sharesA=floor(preValue/priceA(startT,1)/perA)*perA;
    marketA=sharesA*priceA(startT,1);
    tradeFeeA=marketA*(feeA+stampA);
    shortFeeA=marketA*sRateA/cDay;
    ableValue=totalAsset(startT-1,1)-marketH-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA;
    ableRate=ableValue*r/cDay;
    totalAsset(startT,1)=totalAsset(startT-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate;
    accountDetail(startT,:)=[totalAsset(startT,1),marketH,sharesH*tradeFlagH,tradeFeeH,shortFeeH,marketA,sharesA*tradeFlagA,tradeFeeA,shortFeeA,ableValue,ableRate,0,marketA];
    tradeDetail(startT,:)=[sharesH*tradeFlagH,priceH(startT,1),sharesA*tradeFlagA,priceA(startT,1),0,0];
else
    %融券卖开H股，买开A股
    tradeFlagH=-1;
    tradeFlagA=1;
    preValue=-position(startT)*totalAsset(startT-1,1)/(1+feeA+feeH+stampH+sRateH*shortMonth/12);%预留手续费(买开H，卖开A)+融券费用
    sharesA=floor(preValue/priceA(startT,1)/perA)*perA;
    marketA=sharesA*priceA(startT,1);
    tradeFeeA=marketA*feeA;
    shortFeeA=0;
    sharesH=floor(preValue/priceH(startT,1)/perH)*perH;
    marketH=sharesH*priceH(startT,1);
    tradeFeeH=marketH*(feeH+stampH);
    shortFeeH=marketH*sRateH/cDay;%每日计算，交易计提
    ableValue=totalAsset(startT-1,1)-marketA-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA;
    ableRate=ableValue*r/cDay;
    totalAsset(startT,1)=totalAsset(startT-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate;
    accountDetail(startT,:)=[totalAsset(startT,1),marketH,sharesH*tradeFlagH,tradeFeeH,shortFeeH,marketA,sharesA*tradeFlagA,tradeFeeA,shortFeeA,ableValue,ableRate,marketH,0];
    tradeDetail(startT,:)=[sharesH*tradeFlagH,priceH(startT,1),sharesA*tradeFlagA,priceA(startT,1),0,0];
end

% 每天交易模拟
for t=startT+1:totalT
    if totalAsset(t-1,1)<=0
        accountDetail=accountDetail(1:end-1,:);
        return
    end
    if position(t,1)>=0 && position(t-1,1)>=0
        %当前仓位与上一交易日仓位都>=0
        if position(t,1)==position(t-1,1)
            %仓位不变
            sharesH=accountDetail(t-1,3);
            marketH=sharesH*priceH(t,1);
            tradeFeeH=0;
            shortFeeH=0;
            sharesA=-accountDetail(t-1,7);
            marketA=sharesA*priceA(t,1);
            tradeFeeA=0;
            shortFeeA=(dateS(t,1)-dateS(t-1))/cDay*sRateA*accountDetail(t-1,13);
            ableValue=accountDetail(t-1,10)-shortFeeH-shortFeeA;
            ableRate=accountDetail(t-1,10)*r*(dateS(t,1)-dateS(t-1)-1)/cDay+ableValue*r/cDay;
            gainH=(priceH(t,1)-priceH(t-1,1))*sharesH;
            gainA=-(priceA(t,1)-priceA(t-1,1))*sharesA;
            totalAsset(t,1)=totalAsset(t-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate+gainH+gainA;
            accountDetail(t,1:13)=[totalAsset(t,1),marketH,sharesH,tradeFeeH,shortFeeH,marketA,-sharesA,tradeFeeA,shortFeeA,ableValue,ableRate,0,accountDetail(t-1,13)];
            tradeDetail(t,:)=[0,priceH(t,1),0,priceA(t,1),gainH,gainA];
        else
            %仓位有变化
            preValue=position(t,1)*totalAsset(t-1,1)/(1+feeH+feeA+stampA+sRateA*shortMonth/12);%预留手续费(买开H，卖开A)
            sharesH=floor(preValue/priceH(t,1)/perH)*perH;
            marketH=sharesH*priceH(t,1);
            tradeNumH=sharesH-accountDetail(t-1,3);
            
            %买开卖平操作，卖平考虑印花税
            if tradeNumH>=0
                %买开操作
                tradeFeeH=tradeNumH*priceH(t,1)*feeH;
                gainH=accountDetail(t-1,3)*(priceH(t,1)-priceH(t-1,1));                
            else
                % 卖平
                tradeFeeH=-tradeNumH*priceH(t,1)*(feeH+stampH);
                gainH=-tradeNumH*(priceH(t,1)-priceH(t-1,1))+sharesH*(priceH(t,1)-priceH(t-1,1));
            end
            cashH=-tradeNumH*priceH(t,1);
            shortFeeH=0;
            
            sharesA=floor(preValue/priceA(t,1)/perA)*perA;
            marketA=sharesA*priceA(t,1);
            tradeNumA=-sharesA-accountDetail(t-1,7);            
            %买平卖开操作，卖开考虑印花税+融券利息
            if tradeNumA>0
                % 买平
                tradeFeeA=tradeNumA*priceA(t,1)*feeA;
                gainA=-tradeNumA*(priceA(t,1)-priceA(t-1,1))-sharesA*(priceA(t,1)-priceA(t-1,1));
                sHoldMarketA=accountDetail(t-1,13)*sharesA/(-accountDetail(t-1,7));
                sCloseMarketA=accountDetail(t-1,13)*tradeNumA/(-accountDetail(t-1,7));
                shortFeeA=(dateS(t,1)-dateS(t-1))/cDay*sRateA*sHoldMarketA+(dateS(t,1)-dateS(t-1)-1)/cDay*sRateA*sCloseMarketA;
            else
                % 卖开
                tradeFeeA=-tradeNumA*priceA(t,1)*(feeA+stampA);
                gainA=accountDetail(t-1,7)*(priceA(t,1)-priceA(t-1,1));
                sHoldMarketA=accountDetail(t-1,13)-tradeNumA*priceA(t,1);
                shortFeeA=(dateS(t,1)-dateS(t-1))/cDay*sRateA*accountDetail(t-1,13)-tradeNumA*priceA(t,1)*sRateA/cDay;
            end
            cashA=0;
            % 计算交易后状态
            ableValue=accountDetail(t-1,10)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+cashH+cashA;
            ableRate=accountDetail(t-1,10)*r*(dateS(t,1)-dateS(t-1)-1)/cDay+ableValue*r/cDay;
            totalAsset(t,1)=totalAsset(t-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate+gainH+gainA;
            accountDetail(t,1:13)=[totalAsset(t,1),marketH,sharesH,tradeFeeH,shortFeeH,marketA,-sharesA,tradeFeeA,shortFeeA,ableValue,ableRate,0,sHoldMarketA];
            tradeDetail(t,:)=[tradeNumH,priceH(t,1),tradeNumA,priceA(t,1),gainH,gainA];
        end
        
    elseif position(t,1)<0 && position(t-1,1)<0
        %当前仓位与上一交易日仓位都<0
        if position(t,1)==position(t-1,1)
            %仓位不变
            sharesA=accountDetail(t-1,7);
            marketA=sharesA*priceA(t,1);
            tradeFeeA=0;
            shortFeeA=0;
            sharesH=-accountDetail(t-1,3);
            marketH=sharesH*priceH(t,1);
            tradeFeeH=0;
            shortFeeH=(dateS(t,1)-dateS(t-1))/365*sRateH*accountDetail(t-1,12);
            ableValue=accountDetail(t-1,10)-shortFeeH-shortFeeA;
            ableRate=accountDetail(t-1,10)*r*(dateS(t,1)-dateS(t-1)-1)/cDay+ableValue*r/cDay;
            gainA=(priceA(t,1)-priceA(t-1,1))*sharesA;
            gainH=-(priceH(t,1)-priceH(t-1,1))*sharesH;            
            totalAsset(t,1)=totalAsset(t-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate+gainH+gainA;
            accountDetail(t,1:13)=[totalAsset(t,1),marketH,-sharesH,tradeFeeH,shortFeeH,marketA,sharesA,tradeFeeA,shortFeeA,ableValue,ableRate,accountDetail(t-1,12),0];
            tradeDetail(t,:)=[0,priceH(t,1),0,priceA(t,1),gainH,gainA];
        else
            %仓位有变化
            preValue=-position(t,1)*totalAsset(t-1,1)/(1+feeA+feeH+stampH+sRateH*shortMonth/12);%预留手续费(卖开H，买开A)
            sharesA=floor(preValue/priceA(t,1)/perA)*perA;
            marketA=sharesA*priceA(t,1);
            tradeNumA=sharesA-accountDetail(t-1,7);
            
            %买开卖平操作，卖平考虑印花税
            if tradeNumA>=0
                %买开操作
                tradeFeeA=tradeNumA*priceA(t,1)*feeA;
                gainA=accountDetail(t-1,7)*(priceA(t,1)-priceA(t-1,1));                
            else
                % 卖平
                tradeFeeA=-tradeNumA*priceA(t,1)*(feeA+stampA);
                gainA=-tradeNumA*(priceA(t,1)-priceA(t-1,1))+sharesA*(priceA(t,1)-priceA(t-1,1));
            end
            cashA=-tradeNumA*priceA(t,1);
            shortFeeA=0;
            
            sharesH=floor(preValue/priceH(t,1)/perH)*perH;
            marketH=sharesH*priceH(t,1);
            tradeNumH=-sharesH-accountDetail(t-1,3);            
            %买平卖开操作，卖开考虑印花税+融券利息
            if tradeNumH>0
                % 买平
                tradeFeeH=tradeNumH*priceH(t,1)*feeH;
                gainH=-tradeNumH*(priceH(t,1)-priceH(t-1,1))-sharesH*(priceH(t,1)-priceH(t-1,1));
                sHoldMarketH=accountDetail(t-1,12)*sharesH/(-accountDetail(t-1,3));
                sCloseMarketH=accountDetail(t-1,12)*tradeNumH/(-accountDetail(t-1,3));
                shortFeeH=(dateS(t,1)-dateS(t-1))/cDay*sRateH*sHoldMarketH+(dateS(t,1)-dateS(t-1)-1)/cDay*sRateH*sCloseMarketH;
            else
                % 卖开
                tradeFeeH=-tradeNumH*priceH(t,1)*(feeH+stampH);
                gainH=accountDetail(t-1,3)*(priceH(t,1)-priceH(t-1,1));
                sHoldMarketH=accountDetail(t-1,12)-tradeNumH*priceH(t,1);
                shortFeeH=(dateS(t,1)-dateS(t-1))/cDay*sRateH*accountDetail(t-1,12)-tradeNumH*priceH(t,1)*sRateH/cDay;
            end
            cashH=0;
            % 计算交易后状态
            ableValue=accountDetail(t-1,10)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+cashH+cashA;
            ableRate=accountDetail(t-1,10)*r*(dateS(t,1)-dateS(t-1)-1)/cDay+ableValue*r/cDay;
            totalAsset(t,1)=totalAsset(t-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate+gainH+gainA;
            accountDetail(t,1:13)=[totalAsset(t,1),marketH,-sharesH,tradeFeeH,shortFeeH,marketA,sharesA,tradeFeeA,shortFeeA,ableValue,ableRate,sHoldMarketH,0];
            tradeDetail(t,:)=[tradeNumH,priceH(t,1),tradeNumA,priceA(t,1),gainH,gainA];
        end
        
    elseif position(t,1)<0 && position(t-1,1)>=0
        %当前仓位<0,上交易日仓位>=0,即需要先平仓再开仓
        preValue=-position(t,1)*totalAsset(t-1,1)/(1+feeA+feeH+stampH+sRateH*shortMonth/12);%预留手续费(卖开H，买开A)
        
        % A:买平-accountDetail(t-1,7)+买开sharesA
        sharesA=floor(preValue/priceA(t,1)/perA)*perA;
        marketA=sharesA*priceA(t,1);
        tradeNumA=sharesA-accountDetail(t-1,7);        
        tradeFeeA=tradeNumA*priceA(t,1)*feeA;
        gainA=accountDetail(t-1,7)*(priceA(t,1)-priceA(t-1,1));
        cashA=-sharesA*priceA(t,1);
        shortFeeA=(dateS(t,1)-dateS(t-1)-1)/cDay*sRateA*accountDetail(t-1,13);        
        
        %H:卖平accountDetail(t-1,3)+卖开sharesH
        sharesH=floor(preValue/priceH(t,1)/perH)*perH;        
        marketH=sharesH*priceH(t,1);
        tradeNumH=-sharesH-accountDetail(t-1,3);        
        tradeFeeH=-tradeNumH*priceH(t,1)*(feeH+stampH);
        gainH=accountDetail(t-1,3)*(priceH(t,1)-priceH(t-1,1));
        shortFeeH=marketH*sRateH/cDay;
        cashH=accountDetail(t-1,3)*priceH(t,1);
        
        % 计算交易后状态
        ableValue=accountDetail(t-1,10)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+cashH+cashA;
        ableRate=accountDetail(t-1,10)*r*(dateS(t,1)-dateS(t-1)-1)/cDay+ableValue*r/cDay;
        totalAsset(t,1)=totalAsset(t-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate+gainH+gainA;
        accountDetail(t,1:13)=[totalAsset(t,1),marketH,-sharesH,tradeFeeH,shortFeeH,marketA,sharesA,tradeFeeA,shortFeeA,ableValue,ableRate,marketH,0];
        tradeDetail(t,:)=[tradeNumH,priceH(t,1),tradeNumA,priceA(t,1),gainH,gainA];        
        
    elseif position(t,1)>=0 && position(t-1,1)<0
        %当前仓位>=0,上交易日仓位<0,即需要先平仓再开仓
        preValue=position(t,1)*totalAsset(t-1,1)/(1+feeH+feeA+stampA+sRateA*shortMonth/12);%预留手续费(买开H，卖开A)
        
        % H:买平-accountDetail(t-1,3)+买开sharesH
        sharesH=floor(preValue/priceH(t,1)/perH)*perH;
        marketH=sharesH*priceH(t,1);
        tradeNumH=sharesH-accountDetail(t-1,3);        
        tradeFeeH=tradeNumH*priceH(t,1)*feeH;
        gainH=accountDetail(t-1,3)*(priceH(t,1)-priceH(t-1,1));
        cashH=-sharesH*priceH(t,1);
        shortFeeH=(dateS(t,1)-dateS(t-1)-1)/cDay*sRateH*accountDetail(t-1,12);        
        
        %A:卖平accountDetail(t-1,7)+卖开sharesA
        sharesA=floor(preValue/priceA(t,1)/perA)*perA;        
        marketA=sharesA*priceA(t,1);
        tradeNumA=-sharesA-accountDetail(t-1,7);        
        tradeFeeA=-tradeNumA*priceA(t,1)*(feeA+stampA);
        gainA=accountDetail(t-1,7)*(priceA(t,1)-priceA(t-1,1));
        shortFeeA=marketA*sRateA/cDay;
        cashA=accountDetail(t-1,7)*priceA(t,1);
        
        % 计算交易后状态
        ableValue=accountDetail(t-1,10)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+cashH+cashA;
        ableRate=accountDetail(t-1,10)*r*(dateS(t,1)-dateS(t-1)-1)/cDay+ableValue*r/cDay;
        totalAsset(t,1)=totalAsset(t-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate+gainH+gainA;
        accountDetail(t,1:13)=[totalAsset(t,1),marketH,sharesH,tradeFeeH,shortFeeH,marketA,-sharesA,tradeFeeA,shortFeeA,ableValue,ableRate,0,marketA];
        tradeDetail(t,:)=[tradeNumH,priceH(t,1),tradeNumA,priceA(t,1),gainH,gainA];
    end
    
    
end
end