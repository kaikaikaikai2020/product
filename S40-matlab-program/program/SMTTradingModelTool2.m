%% ��װԭ�����򣬱��ں����������ݿ�
%�޸��»�ͼ����
%������������Ȩ��
%ÿ��������ĸ�����ΪȨ�أ�����������0.6�Ŀ������ǣ����Ǿ���0.6��ΪȨ�ء�����
%�����һ����Ʊ0.55�ĸ����ǣ�Ȩ�ؾ�Ϊ0.55.
%����ԭ����ܣ�ԭ�����������һ��
%% SMTTradingModel
% SMTʱ����β���
% ��������--����������
% ���ݴ���--�޳���ָ�ڻ���ǰ��15��������
% ����ģ��--ѡȡ�����ǰM��ƥ��
function [tradeYield,result,tradeDetail,yearDetail,h,assure_ratio] =SMTTradingModelTool2(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate,P)
    %% ��������
    feeOpen=P.feeOpen;
    feeClose=P.feeClose;
    matchRecord=P.matchRecord;%ƥ������Դ������300
    tradeRecord=P.tradeRecord;%��������Դ����ָ�ڻ�������Լ
    tradeMin=P.tradeMin;%ʹ������120����K�����ݽ��з���ƥ��
    dayMin=P.dayMin;%ÿ�������չ�240��1����K��
    M=P.M;%��M����Ϊ���ƵĽ�����
    muchPara=P.muchPara;%�������ǻ��µ�����
    deanMethod=P.deanMethod;%1���ϵ��/2ŷʽ����/3���Ͼ���/4�����پ���
    stopMethod=P.stopMethod;%1���̼�ֹ��/2������ֹ�𣬲������ܷ���/3�������̼�ֹ�𣬷��򴥷���ֹ��
    testStart=P.testStart;
    trade_mode = P.trade_mode;
    if isfield(P,'out_sel')
        out_sel = P.out_sel;
    else
        out_sel = true;
    end
    %endDate='2018-6-1';%���������һ�������գ�Ŀǰ�������µ�2018-5-31
    %% ���ݴ���
    % ���ƥ������Դ�ͽ�������Դ
    [matchOpen,matchClose,matchDate,tradeOpen,tradeClose,tradeDate]=dataSourse(matchRecord,tradeRecord,hsOpen,hsClose,hsDate,ifOpen,ifClose,ifDate);
    %% ����ģ��
    % ƥ�����ݴ���׼��
    outMatchIndex=find(floor(matchDate)>=max(datenum(testStart),floor(ifDate(1))));
    outMatchClose=matchClose(outMatchIndex(1):end);
    outMatchOpen=matchOpen(outMatchIndex(1):end);
    outMatchDate=matchDate(outMatchIndex(1):end);
    
    inMatchDate=matchDate(1:outMatchIndex(1)-1);
    allDay=[1;find(floor(matchDate(1:end-1))~=floor(matchDate(2:end)))+1];
    outDay=[1;find(floor(outMatchDate(1:end-1))~=floor(outMatchDate(2:end)))+1];
    inDay=[1;find(floor(inMatchDate(1:end-1))~=floor(inMatchDate(2:end)))+1];
    % �������ݴ���׼��
    outTradeIndex=find(floor(tradeDate)>=datenum(testStart));
    outTradeClose=tradeClose(outTradeIndex(1):end);
    outTradeOpen=tradeOpen(outTradeIndex(1):end);
    outTradeDate=tradeDate(outTradeIndex(1):end);

    % �����ƥ��
    corrRecord=corrSimilarityMatch2(allDay,outDay,inDay,matchClose,outMatchClose,matchDate,outMatchDate,tradeMin,M,deanMethod);
    % ��ֹ̬����
    [fractalPostion,assure_ratio]=dynamicTrailingStop2(corrRecord,outDay,outMatchDate,outMatchClose,tradeMin,muchPara);

    % �������
    if stopMethod==1
        % 1������ֹ��۵ĵ�һ��K�ߵ����̼�
        tradeYield=tradeSimulate1(fractalPostion,outTradeClose,outTradeDate,feeOpen,feeClose);
    elseif stopMethod==2
        % 2������ֹ����ƽ�֣������Ǽ۸��Ƿ��ܹ�����
        tradeYield=tradeSimulate2(fractalPostion,outTradeClose,outTradeDate,feeOpen,feeClose,dayMin);
    elseif stopMethod==3
        % 3�����ƽ��K�ߵĿ����ѵ���ֹ��ۣ���ʹ�ÿ��̼�ֹ�𣬷���ʹ�õ���K�߸յ���ֹ��۵ĵ�λֹ��
        tradeYield=tradeSimulate3(fractalPostion,outTradeClose,outTradeDate,feeOpen,feeClose,dayMin,outTradeOpen);
    end
    if out_sel
        % ������
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
%% ����˵��
% ƥ����ʷ���ݵ�����ԣ�ɸѡ�����ϵ������һ��ֵ�����ڣ�����¼
% corrRecord{1,1}:����/���ϵ��/��ղ�(1Ϊ��,-1����)/β����������ǵ���
%%
corrRecord =cell(length(outDay),1);
for i=1:length(outDay)
    %disp(['��ʼƥ���������', num2str(i),'�����ݡ���' datestr(floor(outMatchDate(outDay(i))),'yyyy-mm-dd')]);
    corrAll=zeros(length(inDay)+i-1,1);
    for j=1:length(inDay)+i-1
        x=outMatchClose(outDay(i):outDay(i)+tradeMin-1,1)/outMatchClose(outDay(i),1);
        y=matchClose(allDay(j):allDay(j)+tradeMin-1,1)/matchClose(allDay(j),1);
        if deanMethod==1
            % 1���ϵ��
            corr=corrcoef(x,y);
            corrAll(j,1)=corr(1,2);
        elseif deanMethod==2
            % 2ŷ�Ͼ���
            euclidean=sum((x-y).^2);
            corrAll(j,1)=euclidean;
        elseif deanMethod==3            
            % 3 ���Ͼ���
            lsdean=sum(abs(x-y)./(x+y));
            corrAll(j,1)=lsdean;
        elseif deanMethod==4          
            % 4�����پ���
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
        % �жϵ���β�����
        if matchClose(allDay(j)+tradeMin-1,1)-matchClose(allDay(j+1)-1,1)>=0
            % �������¼����Ƿ�
            corrRecord{i,1}(ind_i,3)=-1;
            corrRecord{i,1}(ind_i,4)=max(matchClose(allDay(j)+tradeMin-1:allDay(j+1)-1,1))/matchClose(allDay(j)+tradeMin-1,1)-1;
        elseif matchClose(allDay(j)+tradeMin-1,1)-matchClose(allDay(j+1)-1,1)<0
            % �������¼������
            corrRecord{i,1}(ind_i,3)=1;
            corrRecord{i,1}(ind_i,4)=1-min(matchClose(allDay(j)+tradeMin-1:allDay(j+1)-1,1))/matchClose(allDay(j)+tradeMin-1,1);
        end
    end
end
end
function [matchOpen,matchClose,matchDate,tradeOpen,tradeClose,tradeDate]=dataSourse(matchRecord,tradeRecord,hsOpen,hsClose,hsDate,ifOpen,ifClose,ifDate)
%% ����˵��
% ��������������жϻ���300�͹�ָ�ڻ���ƥ������Դ�ͽ�������Դ
% ���������ݣ�һ�������ڼ���ƥ��ȣ�һ�������ڻز⣨һ���ַ��źţ�һ����ִ�н��ף������ֿ�����ͬ��Ҳ���Բ�ͬ��
%%
% �ж�ƥ������Դ
if matchRecord==1
    matchOpen=hsOpen;
    matchClose=hsClose;
    matchDate=hsDate;
else
    matchOpen=ifOpen;
    matchClose=ifClose;
    matchDate=ifDate;
end
% �жϽ�������Դ
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
%% ����˵��
% ��ʱ�����б�׼��
% ������ֵ�͵�ʱ�����У����Ϊʱ������
%% 
dateVec=datevec(dateSeries);
dateVec(:,5)=dateVec(:,5)+round(dateVec(:,6)/60);
dateVec(:,6)=0;

end

function [fractalPostion,assure_ratio]=dynamicTrailingStop2(corrRecord,outDay,outMatchDate,outMatchClose,tradeMin,muchPara)
%% ����˵��
% ��̬׷��ֹ��:�жϵ�ǰ��������Ƿ񳬹���ʷƥ����ε���������Ⱦ�ֵ���糬������Ϊ���γ���ƫ��
% ���fractalPostion:1��������/2���׷���/3�����±�/4��ʷ��������Ⱦ�ֵ/5����������/6ƽ���±�
%assure_ratio ���նȣ����ǻ��µ��ĸ���
ind_j=1;
fractalPostion(ind_j,1)=0;
assure_ratio = ones(length(outDay)-1,1);
for i=1:length(outDay)
    %disp(['��ʼ�����������', num2str(i),'��Ĳ�λ�������' datestr(floor(outMatchDate(outDay(i))),'yyyy-mm-dd')]);
    lenCorr=size(corrRecord{i,1},1); 
    longIndex=find(corrRecord{i,1}(:,3)==1);
    shortIndex=find(corrRecord{i,1}(:,3)==-1);
    if length(longIndex)/lenCorr>muchPara
        % ��ʷƥ�����--����
        fractalPostion(ind_j,1)=floor(outMatchDate(outDay(i)));%��������
        fractalPostion(ind_j,2)=1;%���׷���
        fractalPostion(ind_j,3)=outDay(i)+tradeMin-1;%�����±�
        lossMeanRate=mean(corrRecord{i,1}(longIndex,4));            
        fractalPostion(ind_j,4)=lossMeanRate;%��ʷ���ۿ�����Ⱦ�ֵ
        if i < length(outDay)
            lossNowRate=1-outMatchClose(outDay(i)+tradeMin:outDay(i+1)-1)./outMatchClose(outDay(i)+tradeMin-1);
        else
            lossNowRate=1-outMatchClose(outDay(i)+tradeMin:end)./outMatchClose(outDay(i)+tradeMin-1);
        end
        indC=find(lossNowRate>=lossMeanRate);
        fractalPostion(ind_j,5)=max(lossNowRate);%�������������
        if isempty(indC)==1
            if i < length(outDay)
                fractalPostion(ind_j,6)=outDay(i+1)-1;%ƽ���±�
            else
                fractalPostion(ind_j,6)=length(outMatchClose);%ƽ���±�
            end
        else
            fractalPostion(ind_j,6)=indC(1)+outDay(i)+tradeMin-1;%ƽ���±�
        end
        assure_ratio(ind_j,1) = length(longIndex)/(length(longIndex)+length(shortIndex));
        ind_j=ind_j+1;
    else
        % ��ʷƥ�����--����
        fractalPostion(ind_j,1)=floor(outMatchDate(outDay(i)));%��������
        fractalPostion(ind_j,2)=-1;%���׷���
        fractalPostion(ind_j,3)=outDay(i)+tradeMin-1;%�����±�
        lossMeanRate=mean(corrRecord{i,1}(shortIndex,4));            
        fractalPostion(ind_j,4)=lossMeanRate;%��ʷ���ۿ�����Ⱦ�ֵ
        if i < length(outDay)
            lossNowRate=outMatchClose(outDay(i)+tradeMin:outDay(i+1)-1)./outMatchClose(outDay(i)+tradeMin-1)-1;
        else
            lossNowRate=outMatchClose(outDay(i)+tradeMin:end)./outMatchClose(outDay(i)+tradeMin-1)-1;
        end
        indC=find(lossNowRate>=lossMeanRate);
        fractalPostion(ind_j,5)=max(lossNowRate);%�������������
        if isempty(indC)==1
            if i < length(outDay)
                fractalPostion(ind_j,6)=outDay(i+1)-1;%ƽ���±�
            else
                fractalPostion(ind_j,6)=length(outMatchClose);%ƽ���±�
            end
        else
            fractalPostion(ind_j,6)=indC(1)+outDay(i)+tradeMin-1;%ƽ���±�
        end
        assure_ratio(ind_j,1) = length(shortIndex)/(length(longIndex)+length(shortIndex));
        ind_j=ind_j+1;
    end
end
end

function [tradeDetail,yearDetail]=outputSMT(tradeYield,outTradeDate)
%% ����˵��
% ���ÿ�ʽ����������ʷ������浽Excel
len=size(tradeYield,1);
tradeDetail{len,1}=[];
for i=1:len
    tradeDetail{i,1}=datestr(tradeYield(i,1),'yyyy-mm-dd');
    tradeDetail{i,2}=datestr(outTradeDate(tradeYield(i,6)),'HH:MM');
    tradeDetail{i,5}=datestr(outTradeDate(tradeYield(i,7)),'HH:MM');
    tradeDetail{i,3}=tradeYield(i,4);
    tradeDetail{i,6}=tradeYield(i,5);    
    if tradeYield(i,3)==1
        tradeDetail{i,4}='��';
    else
        tradeDetail{i,4}='��';
    end
    tradeDetail{i,7}=[num2str(floor(tradeYield(i,2)*10000)/100),'%'];
end
tradeDetail=[{'��������','����ʱ��','���ֵ�λ','���ַ���','ƽ��ʱ��','ƽ�ֵ�λ','������'};tradeDetail];
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
yearDetail=[{'���','���״���','�ۼ�������','�ɹ���','���س�','ӯ������','�������'};yearDetail];
end
function [result1,h]=outResultCaculate(tradeYield,feeOpen,feeClose,trade_mode)
%% ����˵��
% �����ǰ���׽����ͳ������
%%
allDayNum=length(tradeYield(:,1));
longNum=length(find(tradeYield(:,3)==1));
shortNum=length(find(tradeYield(:,3)==-1));
tradeNum=longNum+shortNum;
winNum=length(find(tradeYield(:,2)>0));
loseNum=length(find(tradeYield(:,2)<0));
winRate=mean(tradeYield(tradeYield(:,2)>0,2));
loseRate=mean(tradeYield(tradeYield(:,2)<0,2));
%ֻ��� ����
if eq(trade_mode,1)
    Yield=tradeYield(:,2);
    Yield(tradeYield(:,3)<0)=0;
    % �����ǽ��׷��õ��ۼ�������
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

result_info = {'�ܼƽ�����','���״���','ӯ������','�������','ƽ��ӯ����%','ƽ��������%',...
    '�ɹ���%','�ۼ�������%','���س���%','�껯������%','�ۼ����棨�������ѣ�%'};
result(5:end) = result(5:end)*100;
result1 = [result_info;num2cell(result)]';

%% ������
% disp('���������£�')
% disp(['�ܼƽ����գ�' num2str(allDayNum)])
% disp(['���״�����' num2str(tradeNum)])
% disp(['ӯ��������' num2str(winNum)])
% disp(['���������' num2str(loseNum)])
% disp(['ƽ��ӯ���ʣ�' num2str(floor(winRate*10000)/100) '%'])
% disp(['ƽ�������ʣ�' num2str(floor(loseRate*10000)/100) '%'])
% disp(['�ɹ��ʣ�' num2str(floor(winNum/tradeNum*10000)/100) '%'])
% disp(['�ۼ������ʣ�' num2str(floor(cumYield(end)*10000)/100) '%'])
% disp(['���س��ʣ�' num2str(floor(maxDrawDown*10000)/100) '%'])
% disp(['�껯�����ʣ�' num2str(floor(yearYield*10000)/100) '%'])

%% ��ͼ
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
legend('�н��׷�','�޽��׷�','location','best')
title({['�ۼ�����' num2str(floor(cumYield(end)*10000)/100) '%���껯����' num2str(floor(yearYield*10000)/100) '%��' '���س�' num2str(floor(maxDrawDown*10000)/100) '%' ],...
    ['�ܽ��״���' num2str(tradeNum) '��ƽ��ӯ����' num2str(floor(winRate*10000)/100) '%��' 'ƽ��������' num2str(floor(loseRate*10000)/100) '%']})

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
title('�н��׷�������������')

end
function tradeYield=tradeSimulate1(fractalPostion,outTradeClose,outTradeDate,feeOpen,feeClose)
%% ����˵��
% ���ݸ����Ĳ�λ�Թ�ָ���н���
% ֹ���������ʹ�õ���ֹ��۵ĵ�һ��K�ߵ����̼۽���ֹ�����
% ������tradeYield:����/��������/���տ��ַ���/���ֵ�λ/ƽ�ֵ�λ/����ʱ��/ƽ��ʱ��
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
%% ����˵��
% ���ݸ����Ĳ�λ�Թ�ָ���н���
% ֹ�������������ֹ����ƽ�֣������Ǽ۸��Ƿ��ܹ�����
% ������tradeYield:����/��������/���տ��ַ���/���ֵ�λ/ƽ�ֵ�λ/����ʱ��/ƽ��ʱ��
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
%% ����˵��
% ���ݸ����Ĳ�λ�Թ�ָ���н���
% ֹ������������ƽ��K�ߵĿ����ѵ���ֹ��ۣ���ʹ�ÿ��̼�ֹ�𣬷���ʹ�õ���K�߸յ���ֹ��۵ĵ�λֹ��
% ������tradeYield:����/��������/���տ��ַ���/���ֵ�λ/ƽ�ֵ�λ/����ʱ��/ƽ��ʱ��
%%
outTradeDay=unique(floor(outTradeDate));
tradeNum=length(outTradeDay);
tradeYield(1,1)=0;
ind=find(outTradeDay==fractalPostion(1,1));
%������ͳ�Ʒ�������
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
