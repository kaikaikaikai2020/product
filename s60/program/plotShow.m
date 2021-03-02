function plotShow(priceH,priceA,Xt,position,accountDetail,dateS,stockName,result)
%% 函数说明
% 直观展示股票对及价差的走势
%%
figure
lenT=length(accountDetail(:,1));
dateS=dateS(1:lenT);
priceH=priceH(1:lenT);
priceA=priceA(1:lenT);
Xt=Xt(1:lenT);
position=position(1:lenT);
subplot(2,2,1)
plot(priceH,'r')
hold on
plot(priceA,'b')
legend('H股票价格','A股票价格','location','best')
set(gca,'XTick',1:floor(lenT/5):lenT);
set(gca,'XTickLabel',datestr(dateS(1:floor(lenT/5):lenT),'yyyymmdd'));
title('H股与A股价格走势')
subplot(2,2,2)
plot(accountDetail(:,1),'r')
title(['累计收益率' num2str(floor(result(1,1)*100)) '%,年化收益' num2str(floor(result(1,2)*100)) '%,最大回撤'  num2str(floor(result(1,8)*100)) '%,夏普比率'  num2str(floor(result(1,4)*100)/100)])
set(gca,'XTick',1:floor(lenT/5):lenT);
set(gca,'XTickLabel',datestr(dateS(1:floor(lenT/5):lenT),'yyyymmdd'));
subplot(2,2,3)
bar(position,'b')
hold on
plot(Xt,'r')
legend('持仓仓位','lnH-lnA价差走势','location','best')
title('仓位与价差')
set(gca,'XTick',1:floor(lenT/5):lenT);
set(gca,'XTickLabel',datestr(dateS(1:floor(lenT/5):lenT),'yyyymmdd'));
subplot(2,2,4)
dateV=datevec(dateS);
index=[1;find(diff(dateV(:,1)));length(dateS)];
yearS=dateV(index(2:end));
yearRate=accountDetail(index(2:end),1)./accountDetail(index(1:end-1),1)-1;
redRate=yearRate;
redRate(redRate<=0)=0;
blueRate=yearRate;
blueRate(blueRate>0)=0;
bar(redRate,'r')
hold on
bar(blueRate,'b')
set(gca,'xticklabel',yearS)
title('年度收益情况')
%suptitle(['"',stockName,'"','配对交易结果'])
%title(['"',stockName,'"','配对交易结果'])
end