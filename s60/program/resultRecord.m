function result=resultRecord(accountDetail,tDay)
%% ����˵��
% ������Խ��ײ���¼���
% result:1������/2�껯������/3�겨����/4�����ձ���/5ƽ����������/6�ձ�׼��/7�����ձ���/8���س�/9������ռ��/10��ȯ����ռ��
%%
Rf=3/100;%�������ձ��ʵ��޷�������
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