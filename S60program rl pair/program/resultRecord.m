function result=resultRecord(accountDetail,tDay)
%% 函数说明
% 分析配对交易并记录结果
% result:1收益率/2年化收益率/3年波动率/4年夏普比率/5平均日收益率/6日标准差/7日夏普比率/8最大回撤/9手续费占比/10融券费用占比
%%
Rf=3/100;%计算夏普比率的无风险利率
yieldS=accountDetail(:,1);
totalT=length(yieldS);
rateS=[0;yieldS(2:end,1)./yieldS(1:end-1,1)-1];
winRate=yieldS(end,1)/yieldS(1,1)-1;
yearRate=winRate*tDay/totalT;
try
maxDrawDown=maxdrawdown(yieldS);
catch
    maxDrawDown=100;
end
yStd=std(rateS)*sqrt(tDay);
sharpRatio=(yearRate-Rf)/yStd;
meanRate=mean(rateS);
dStd=std(rateS);
dSharpRatio=(meanRate-Rf/tDay)/dStd;
feeS=cumsum(sum(accountDetail(:,[4,8]),2));
shortS=cumsum(sum(accountDetail(:,[5,9]),2));
feeRatio=feeS(end)/yieldS(1,1);
shortRatio=shortS(end)/yieldS(1,1);
result=[winRate,yearRate,yStd,sharpRatio,meanRate,dStd,dSharpRatio,maxDrawDown,feeRatio,shortRatio];
end