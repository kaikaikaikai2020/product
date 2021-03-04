function [accountDetail,tradeDetail]=tradeSimulate(position,priceH,priceA,dateS,V0,feeH,feeA,sRateH,sRateA,stampH,stampA,r,cDay)
%% ����˵��
% ���ݼ�����������ò�λ��ÿ�նԹ�Ʊ�Խ���ģ�⽻�ף�����ģ����ʵ����
% ��ʼͶ�ʽ��ΪV0,��ȯ������ͷ����β������Ȼ��������
% �޷�����Ϣ��ͷ����β��������Ȼ����������
% accountDetail:1���ʲ�/2H��ֵ/3H�ֲ�����������/4H���׷���/5H��ȯ��Ϣ/6A��ֵ/7A�ֲ�����������/8A���׷���/9A��ȯ��Ϣ/10�����ʽ�/11�޷�����Ϣ/12H��ȯ���/13A��ȯ���
% tradeDetail:1H��������������/2H���׼۸�/3A��������������/4A���׼۸�/5H����ӯ��/6A����ӯ��
%%
perH=100;
perA=100;
shortMonth=6;%Ԥ��6������ȯ���ã���ȯ�����Ϊ6���£�
totalT=length(dateS);
startT=min(find(position~=0));
totalAsset(1:startT-1,1)=V0;
accountDetail(1:startT-1,13)=0;
accountDetail(1:startT-1,[1,10])=V0;%����δ����ǰ�����޷�����Ϣ
tradeDetail(1:startT-1,1:6)=0;

% ��ʼ����
if position(startT)>0
    %��H�ɣ���ȯ����A��
    tradeFlagH=1;
    tradeFlagA=-1;
    preValue=position(startT)*totalAsset(startT-1,1)/(1+feeH+feeA+stampA+sRateA*shortMonth/12);%Ԥ��������(��H������A)+shortMonth���µ���ȯ����
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
    %��ȯ����H�ɣ���A��
    tradeFlagH=-1;
    tradeFlagA=1;
    preValue=-position(startT)*totalAsset(startT-1,1)/(1+feeA+feeH+stampH+sRateH*shortMonth/12);%Ԥ��������(��H������A)+��ȯ����
    sharesA=floor(preValue/priceA(startT,1)/perA)*perA;
    marketA=sharesA*priceA(startT,1);
    tradeFeeA=marketA*feeA;
    shortFeeA=0;
    sharesH=floor(preValue/priceH(startT,1)/perH)*perH;
    marketH=sharesH*priceH(startT,1);
    tradeFeeH=marketH*(feeH+stampH);
    shortFeeH=marketH*sRateH/cDay;%ÿ�ռ��㣬���׼���
    ableValue=totalAsset(startT-1,1)-marketA-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA;
    ableRate=ableValue*r/cDay;
    totalAsset(startT,1)=totalAsset(startT-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate;
    accountDetail(startT,:)=[totalAsset(startT,1),marketH,sharesH*tradeFlagH,tradeFeeH,shortFeeH,marketA,sharesA*tradeFlagA,tradeFeeA,shortFeeA,ableValue,ableRate,marketH,0];
    tradeDetail(startT,:)=[sharesH*tradeFlagH,priceH(startT,1),sharesA*tradeFlagA,priceA(startT,1),0,0];
end

% ÿ�콻��ģ��
for t=startT+1:totalT
    if totalAsset(t-1,1)<=0
        accountDetail=accountDetail(1:end-1,:);
        return
    end
    if position(t,1)>=0 && position(t-1,1)>=0
        %��ǰ��λ����һ�����ղ�λ��>=0
        if position(t,1)==position(t-1,1)
            %��λ����
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
            %��λ�б仯
            preValue=position(t,1)*totalAsset(t-1,1)/(1+feeH+feeA+stampA+sRateA*shortMonth/12);%Ԥ��������(��H������A)
            sharesH=floor(preValue/priceH(t,1)/perH)*perH;
            marketH=sharesH*priceH(t,1);
            tradeNumH=sharesH-accountDetail(t-1,3);
            
            %����ƽ��������ƽ����ӡ��˰
            if tradeNumH>=0
                %�򿪲���
                tradeFeeH=tradeNumH*priceH(t,1)*feeH;
                gainH=accountDetail(t-1,3)*(priceH(t,1)-priceH(t-1,1));                
            else
                % ��ƽ
                tradeFeeH=-tradeNumH*priceH(t,1)*(feeH+stampH);
                gainH=-tradeNumH*(priceH(t,1)-priceH(t-1,1))+sharesH*(priceH(t,1)-priceH(t-1,1));
            end
            cashH=-tradeNumH*priceH(t,1);
            shortFeeH=0;
            
            sharesA=floor(preValue/priceA(t,1)/perA)*perA;
            marketA=sharesA*priceA(t,1);
            tradeNumA=-sharesA-accountDetail(t-1,7);            
            %��ƽ������������������ӡ��˰+��ȯ��Ϣ
            if tradeNumA>0
                % ��ƽ
                tradeFeeA=tradeNumA*priceA(t,1)*feeA;
                gainA=-tradeNumA*(priceA(t,1)-priceA(t-1,1))-sharesA*(priceA(t,1)-priceA(t-1,1));
                sHoldMarketA=accountDetail(t-1,13)*sharesA/(-accountDetail(t-1,7));
                sCloseMarketA=accountDetail(t-1,13)*tradeNumA/(-accountDetail(t-1,7));
                shortFeeA=(dateS(t,1)-dateS(t-1))/cDay*sRateA*sHoldMarketA+(dateS(t,1)-dateS(t-1)-1)/cDay*sRateA*sCloseMarketA;
            else
                % ����
                tradeFeeA=-tradeNumA*priceA(t,1)*(feeA+stampA);
                gainA=accountDetail(t-1,7)*(priceA(t,1)-priceA(t-1,1));
                sHoldMarketA=accountDetail(t-1,13)-tradeNumA*priceA(t,1);
                shortFeeA=(dateS(t,1)-dateS(t-1))/cDay*sRateA*accountDetail(t-1,13)-tradeNumA*priceA(t,1)*sRateA/cDay;
            end
            cashA=0;
            % ���㽻�׺�״̬
            ableValue=accountDetail(t-1,10)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+cashH+cashA;
            ableRate=accountDetail(t-1,10)*r*(dateS(t,1)-dateS(t-1)-1)/cDay+ableValue*r/cDay;
            totalAsset(t,1)=totalAsset(t-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate+gainH+gainA;
            accountDetail(t,1:13)=[totalAsset(t,1),marketH,sharesH,tradeFeeH,shortFeeH,marketA,-sharesA,tradeFeeA,shortFeeA,ableValue,ableRate,0,sHoldMarketA];
            tradeDetail(t,:)=[tradeNumH,priceH(t,1),tradeNumA,priceA(t,1),gainH,gainA];
        end
        
    elseif position(t,1)<0 && position(t-1,1)<0
        %��ǰ��λ����һ�����ղ�λ��<0
        if position(t,1)==position(t-1,1)
            %��λ����
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
            %��λ�б仯
            preValue=-position(t,1)*totalAsset(t-1,1)/(1+feeA+feeH+stampH+sRateH*shortMonth/12);%Ԥ��������(����H����A)
            sharesA=floor(preValue/priceA(t,1)/perA)*perA;
            marketA=sharesA*priceA(t,1);
            tradeNumA=sharesA-accountDetail(t-1,7);
            
            %����ƽ��������ƽ����ӡ��˰
            if tradeNumA>=0
                %�򿪲���
                tradeFeeA=tradeNumA*priceA(t,1)*feeA;
                gainA=accountDetail(t-1,7)*(priceA(t,1)-priceA(t-1,1));                
            else
                % ��ƽ
                tradeFeeA=-tradeNumA*priceA(t,1)*(feeA+stampA);
                gainA=-tradeNumA*(priceA(t,1)-priceA(t-1,1))+sharesA*(priceA(t,1)-priceA(t-1,1));
            end
            cashA=-tradeNumA*priceA(t,1);
            shortFeeA=0;
            
            sharesH=floor(preValue/priceH(t,1)/perH)*perH;
            marketH=sharesH*priceH(t,1);
            tradeNumH=-sharesH-accountDetail(t-1,3);            
            %��ƽ������������������ӡ��˰+��ȯ��Ϣ
            if tradeNumH>0
                % ��ƽ
                tradeFeeH=tradeNumH*priceH(t,1)*feeH;
                gainH=-tradeNumH*(priceH(t,1)-priceH(t-1,1))-sharesH*(priceH(t,1)-priceH(t-1,1));
                sHoldMarketH=accountDetail(t-1,12)*sharesH/(-accountDetail(t-1,3));
                sCloseMarketH=accountDetail(t-1,12)*tradeNumH/(-accountDetail(t-1,3));
                shortFeeH=(dateS(t,1)-dateS(t-1))/cDay*sRateH*sHoldMarketH+(dateS(t,1)-dateS(t-1)-1)/cDay*sRateH*sCloseMarketH;
            else
                % ����
                tradeFeeH=-tradeNumH*priceH(t,1)*(feeH+stampH);
                gainH=accountDetail(t-1,3)*(priceH(t,1)-priceH(t-1,1));
                sHoldMarketH=accountDetail(t-1,12)-tradeNumH*priceH(t,1);
                shortFeeH=(dateS(t,1)-dateS(t-1))/cDay*sRateH*accountDetail(t-1,12)-tradeNumH*priceH(t,1)*sRateH/cDay;
            end
            cashH=0;
            % ���㽻�׺�״̬
            ableValue=accountDetail(t-1,10)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+cashH+cashA;
            ableRate=accountDetail(t-1,10)*r*(dateS(t,1)-dateS(t-1)-1)/cDay+ableValue*r/cDay;
            totalAsset(t,1)=totalAsset(t-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate+gainH+gainA;
            accountDetail(t,1:13)=[totalAsset(t,1),marketH,-sharesH,tradeFeeH,shortFeeH,marketA,sharesA,tradeFeeA,shortFeeA,ableValue,ableRate,sHoldMarketH,0];
            tradeDetail(t,:)=[tradeNumH,priceH(t,1),tradeNumA,priceA(t,1),gainH,gainA];
        end
        
    elseif position(t,1)<0 && position(t-1,1)>=0
        %��ǰ��λ<0,�Ͻ����ղ�λ>=0,����Ҫ��ƽ���ٿ���
        preValue=-position(t,1)*totalAsset(t-1,1)/(1+feeA+feeH+stampH+sRateH*shortMonth/12);%Ԥ��������(����H����A)
        
        % A:��ƽ-accountDetail(t-1,7)+��sharesA
        sharesA=floor(preValue/priceA(t,1)/perA)*perA;
        marketA=sharesA*priceA(t,1);
        tradeNumA=sharesA-accountDetail(t-1,7);        
        tradeFeeA=tradeNumA*priceA(t,1)*feeA;
        gainA=accountDetail(t-1,7)*(priceA(t,1)-priceA(t-1,1));
        cashA=-sharesA*priceA(t,1);
        shortFeeA=(dateS(t,1)-dateS(t-1)-1)/cDay*sRateA*accountDetail(t-1,13);        
        
        %H:��ƽaccountDetail(t-1,3)+����sharesH
        sharesH=floor(preValue/priceH(t,1)/perH)*perH;        
        marketH=sharesH*priceH(t,1);
        tradeNumH=-sharesH-accountDetail(t-1,3);        
        tradeFeeH=-tradeNumH*priceH(t,1)*(feeH+stampH);
        gainH=accountDetail(t-1,3)*(priceH(t,1)-priceH(t-1,1));
        shortFeeH=marketH*sRateH/cDay;
        cashH=accountDetail(t-1,3)*priceH(t,1);
        
        % ���㽻�׺�״̬
        ableValue=accountDetail(t-1,10)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+cashH+cashA;
        ableRate=accountDetail(t-1,10)*r*(dateS(t,1)-dateS(t-1)-1)/cDay+ableValue*r/cDay;
        totalAsset(t,1)=totalAsset(t-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate+gainH+gainA;
        accountDetail(t,1:13)=[totalAsset(t,1),marketH,-sharesH,tradeFeeH,shortFeeH,marketA,sharesA,tradeFeeA,shortFeeA,ableValue,ableRate,marketH,0];
        tradeDetail(t,:)=[tradeNumH,priceH(t,1),tradeNumA,priceA(t,1),gainH,gainA];        
        
    elseif position(t,1)>=0 && position(t-1,1)<0
        %��ǰ��λ>=0,�Ͻ����ղ�λ<0,����Ҫ��ƽ���ٿ���
        preValue=position(t,1)*totalAsset(t-1,1)/(1+feeH+feeA+stampA+sRateA*shortMonth/12);%Ԥ��������(��H������A)
        
        % H:��ƽ-accountDetail(t-1,3)+��sharesH
        sharesH=floor(preValue/priceH(t,1)/perH)*perH;
        marketH=sharesH*priceH(t,1);
        tradeNumH=sharesH-accountDetail(t-1,3);        
        tradeFeeH=tradeNumH*priceH(t,1)*feeH;
        gainH=accountDetail(t-1,3)*(priceH(t,1)-priceH(t-1,1));
        cashH=-sharesH*priceH(t,1);
        shortFeeH=(dateS(t,1)-dateS(t-1)-1)/cDay*sRateH*accountDetail(t-1,12);        
        
        %A:��ƽaccountDetail(t-1,7)+����sharesA
        sharesA=floor(preValue/priceA(t,1)/perA)*perA;        
        marketA=sharesA*priceA(t,1);
        tradeNumA=-sharesA-accountDetail(t-1,7);        
        tradeFeeA=-tradeNumA*priceA(t,1)*(feeA+stampA);
        gainA=accountDetail(t-1,7)*(priceA(t,1)-priceA(t-1,1));
        shortFeeA=marketA*sRateA/cDay;
        cashA=accountDetail(t-1,7)*priceA(t,1);
        
        % ���㽻�׺�״̬
        ableValue=accountDetail(t-1,10)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+cashH+cashA;
        ableRate=accountDetail(t-1,10)*r*(dateS(t,1)-dateS(t-1)-1)/cDay+ableValue*r/cDay;
        totalAsset(t,1)=totalAsset(t-1,1)-tradeFeeH-shortFeeH-tradeFeeA-shortFeeA+ableRate+gainH+gainA;
        accountDetail(t,1:13)=[totalAsset(t,1),marketH,sharesH,tradeFeeH,shortFeeH,marketA,-sharesA,tradeFeeA,shortFeeA,ableValue,ableRate,0,marketA];
        tradeDetail(t,:)=[tradeNumH,priceH(t,1),tradeNumA,priceA(t,1),gainH,gainA];
    end
    
    
end
end