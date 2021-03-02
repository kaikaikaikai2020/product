function plotShow(priceH,priceA,Xt,position,accountDetail,dateS,stockName,result)
%% ����˵��
% ֱ��չʾ��Ʊ�Լ��۲������
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
legend('H��Ʊ�۸�','A��Ʊ�۸�','location','best')
set(gca,'XTick',1:floor(lenT/5):lenT);
set(gca,'XTickLabel',datestr(dateS(1:floor(lenT/5):lenT),'yyyymmdd'));
title('H����A�ɼ۸�����')
subplot(2,2,2)
plot(accountDetail(:,1),'r')
title(['�ۼ�������' num2str(floor(result(1,1)*100)) '%,�껯����' num2str(floor(result(1,2)*100)) '%,���س�'  num2str(floor(result(1,8)*100)) '%,���ձ���'  num2str(floor(result(1,4)*100)/100)])
set(gca,'XTick',1:floor(lenT/5):lenT);
set(gca,'XTickLabel',datestr(dateS(1:floor(lenT/5):lenT),'yyyymmdd'));
subplot(2,2,3)
bar(position,'b')
hold on
plot(Xt,'r')
legend('�ֲֲ�λ','lnH-lnA�۲�����','location','best')
title('��λ��۲�')
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
title('����������')
%suptitle(['"',stockName,'"','��Խ��׽��'])
%title(['"',stockName,'"','��Խ��׽��'])
end