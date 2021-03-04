%{
1���۲ⴰ��������һ�������
2������β����close�۸�����
3��ͳ�Ƶڶ����������ȷ������
4���Եڶ���close�۸�����
5����������stop loss����
6�������ѻ���˫��ǧ��֮1.5
%}
%% SMTTradingModel
% SMTʱ����β���
% ��������--����������
% ���ݴ���--�޳���ָ�ڻ���ǰ��15��������
% ����ģ��--ѡȡ�����ǰM��ƥ��
function fractalPostion =SMTTradingModelTool_update1(hsClose,hsOpen,hsDate,ifClose,ifOpen,ifDate,P)
    %% ��������
    matchRecord=P.matchRecord;%ƥ������Դ������300
    tradeRecord=P.tradeRecord;%��������Դ����ָ�ڻ�������Լ
    tradeMin=P.tradeMin;%ʹ������120����K�����ݽ��з���ƥ��
    M=P.M;%��M����Ϊ���ƵĽ�����
    muchPara=P.muchPara;%�������ǻ��µ�����
    deanMethod=P.deanMethod;%1���ϵ��/2ŷʽ����/3���Ͼ���/4�����پ���
    testStart=P.testStart;
    cut_return = P.cut_return;
    %endDate='2018-6-1';%���������һ�������գ�Ŀǰ�������µ�2018-5-31
    %% ���ݴ���
    % ���ƥ������Դ�ͽ�������Դ
    [~,matchClose,matchDate]=dataSourse(matchRecord,tradeRecord,hsOpen,hsClose,hsDate,ifOpen,ifClose,ifDate);
    %% ����ģ��
    % ƥ�����ݴ���׼��
    outMatchIndex=find(floor(matchDate)>=max(datenum(testStart),floor(ifDate(1))));
    outMatchClose=matchClose(outMatchIndex(1):end);
    outMatchDate=matchDate(outMatchIndex(1):end);
    
    inMatchDate=matchDate(1:outMatchIndex(1)-1);
    allDay=[1;find(floor(matchDate(1:end-1))~=floor(matchDate(2:end)))+1];
    outDay=[1;find(floor(outMatchDate(1:end-1))~=floor(outMatchDate(2:end)))+1];
    inDay=[1;find(floor(inMatchDate(1:end-1))~=floor(inMatchDate(2:end)))+1];
    %�����ƥ��
    %corrRecord{1,1}:1 ����/2 ���ϵ��/3 ��ղ�(1Ϊ��,-1����)/ 4 β����������ǵ���  5����������ǵ�
    corrRecord=corrSimilarityMatch2(allDay,outDay,inDay,matchClose,outMatchClose,matchDate,tradeMin,M,deanMethod);
    % �����ź�
    fractalPostion=dynamicTrailingStop2(corrRecord,outDay,outMatchDate,outMatchClose,tradeMin,muchPara,cut_return);
    fractalPostion = fractalPostion(:,1:2);
end

function corrRecord=corrSimilarityMatch2(allDay,outDay,inDay,matchClose,outMatchClose,matchDate,tradeMin,M,deanMethod)
%% ����˵��
% ƥ����ʷ���ݵ�����ԣ�ɸѡ�����ϵ������һ��ֵ�����ڣ�����¼
% corrRecord{1,1}:1 ����/2 ���ϵ��/3 ��ղ�(1Ϊ��,-1����)/ 4 β����������ǵ���  5����������ǵ�
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
        elseif deanMethod==5          
            % 5dynamic time warping ����
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
        % �ж�����β�����
        if matchClose(allDay(j)+tradeMin,1)-matchClose(allDay(j+1)+tradeMin-1,1)>=0
            % ����
            corrRecord{i,1}(ind_i,3)=-1;            
        elseif matchClose(allDay(j)+tradeMin,1)-matchClose(allDay(j+1)+tradeMin-1,1)<0
            % ����
            corrRecord{i,1}(ind_i,3)=1;            
        end
        corrRecord{i,1}(ind_i,5) = matchClose(allDay(j+1)+tradeMin-1,1)/matchClose(allDay(j)+tradeMin,1)-1;
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

function fractalPostion=dynamicTrailingStop2(corrRecord,outDay,outMatchDate,outMatchClose,tradeMin,muchPara,cut_return)
%% ����˵��
% ��̬׷��ֹ��:�жϵ�ǰ��������Ƿ񳬹���ʷƥ����ε���������Ⱦ�ֵ���糬������Ϊ���γ���ƫ��
% ���fractalPostion:1��������/2���׷���/3�����±�/4��ʷ��������Ⱦ�ֵ/5����������/6ƽ���±�
ind_j=1;
fractalPostion(ind_j,1)=0;
for i=1:length(outDay)
    %disp(['��ʼ�����������', num2str(i),'��Ĳ�λ�������' datestr(floor(outMatchDate(outDay(i))),'yyyy-mm-dd')]);
    lenCorr=size(corrRecord{i,1},1); 
    longIndex=find(corrRecord{i,1}(:,3)==1);
    %������һ���ж�
    if length(longIndex)/lenCorr>muchPara
        % ��ʷƥ�����--����
        fractalPostion(ind_j,1)=floor(outMatchDate(outDay(i)));%��������
        fractalPostion(ind_j,2)=1;%���׷���
        fractalPostion(ind_j,3)=outDay(i)+tradeMin-1;%�����±�
        ind_j=ind_j+1;
    else
        % ��ʷƥ�����--����
        fractalPostion(ind_j,1)=floor(outMatchDate(outDay(i)));%��������
        fractalPostion(ind_j,2)=-1;%���׷���
        fractalPostion(ind_j,3)=outDay(i)+tradeMin-1;%�����±�
        ind_j=ind_j+1;
    end
end
end


