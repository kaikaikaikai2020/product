function [v,v_str,sta_val] = curve_static(y,wid_year,disp_sel)
if nargin <2
    wid_year = [];
end
if isempty(wid_year)
    wid_year = 244;
end
if nargin < 3
    disp_sel = [];
end

if isempty(disp_sel)
    disp_sel = false;
end
[m,n] = size(y);
if eq(m,1)
    y = y';
end
% y = cumprod(1+rand(1000,1)/1000);
%(AC3277/100)^(244/COUNT(AC120:AC3277))-1
%1�껯������
v_str{1} = '�껯������';
v(1) = (y(end)/y(1))^(wid_year/length(y))-1;
sta_val.nh = v(1);
%2���س�
v_str{2} = '���س�';
v(2) = min(getdrawdown(y));
sta_val.drawdown = v(2);
%3-5����һ����С\���\ƽ��ƽ������
if length(y)>=wid_year
    v1 = y(wid_year:end)./y(1:end-wid_year+1)-1;
else
    v1 = nan;
end
v_str{3} = '����1���������';
v(3) = min(v1);
v_str{4} = '����1���������';
v(4) = max(v1);
v_str{5} = '����1��ƽ������';
v(5) = mean(v1);
sta_val.min_return = v(3);
sta_val.max_return = v(4);
sta_val.ave_return = v(5);
%6������������� ���������ǲ�������ӯ�����������
testv1 = find(diff(y)>0);
    %�ҳ��������ǵ������
testv2 = diff(testv1);
    %��������Ų�Ϊ1�����ߴ���1
testv3 = [(1:size(testv2))',testv2];
testv3(eq(testv2,1),:)=[];
v_str{6} = '���������������';
if  size(testv3,1)>1
    v(6) = max(diff(testv3(:,1)));
else
    v(6) = 1;
end
sta_val.max_suc_up = v(6);
%7��������µ���
testv1 = find(diff(y)<0);
    %�ҳ��µ������
testv2 = diff(testv1);
    %��������Ų�Ϊ1�����ߴ���1
testv3 = [(1:size(testv2))',testv2];
testv3(eq(testv2,1),:)=[];
v_str{7} = '��������µ�����';
if size(testv3,1)>1
    v(7) = max(diff(testv3(:,1)));
else
    v(7) = 1;
end
sta_val.max_suc_down = v(7);
%8���沨����
    %std(ÿ��������)*sqrt(252)
v_str{8} = '���沨����';
v(8) = std(y(2:end)./y(1:end-1)-1)*sqrt(wid_year);%�껯��׼��
sta_val.std_return = v(8);
%9���ձ��ʣ��������������ʵ�252��ƽ��/�����������ʵ�252�ձ�׼�*[252^��1/2��]
%v(9) = v(1)/v(8)*sqrt(252);
v_str{9} = '���ձ���';
%v(9) = (((mean(y(2:end)./y(1:end-1)-1)-3/100/252))/std(y(2:end)./y(1:end-1)-1))*sqrt(252);
temp = y(2:end)./y(1:end-1)-1;
temp(isinf(temp)|isnan(temp)) = [];
%v(9) = ((mean(temp)-(exp(log(1.03)/252)-1)))/(std(temp))*sqrt(wid_year);
v(9) = (mean(temp))/(std(temp))*sqrt(wid_year); %20200504���� �޷��������Ϊ0
sta_val.sharp = v(9);
%10��¸߼��
v_str{10} = '��¸߼��';
temp = getdrawdown(y);
ind = find(eq(temp,0));
if length(ind)>1
    v(10) = max(diff(ind));
else
    v(10) = nan;
end
sta_val.max_inter = v(10);
%11 С��0������Ĳ�����
v_str{11} = '�µ�������';
v(11) = std(y(2:end)./y(1:end-1)-1);
sta_val.down_std = v(11);
v1 = y(wid_year:end)./y(1:end-wid_year+1)-1;
v_str{12} = '����1��������λ��';
v(12) = median(v1);
sta_val.median_return = v(12);
v_str{13} = '����1������90%��λ��';
v(13) = prctile(v1,90);
sta_val.ninety_quantile = v(13);

v_str{14} = '��Ϣ����';
v(14) = v(1)/v(8);
sta_val.ninety_quantile = v(14);

temp = y(2:end)./y(1:end-1)-1;
temp(eq(temp,0)) = [];
v_str{15} = 'ʤ��';
v(15) = sum(temp>0)/length(temp);
sta_val.ninety_quantile = v(14);


if disp_sel
str = [];
for i = 1:length(v)
    str = [str,sprintf('%s: %0.4f\n',v_str{i},v(i))];
end
sprintf('�ز����߲�����%s \n',str)
end
end


function v = getdrawdown(x)
    v = zeros(size(x));
    for i = 1:size(x,2)
        v(:,i) = x(:,i)./cummax(x(:,i))-1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Info
% %�������۲���
% %�껯����-���㹫ʽΪ�����˻����ռ�ֵ/�˻���ʼ��ֵ��^��250/�ز��ڼ���������-1
% value = [];
% value.nianhuashouyi = (MT(end)/MT(1))^(250/(N*Wnum*7))-1;
% %���س�-���㹫ʽΪ��min(�˻����ռ�ֵ / ����֮ǰ�˻���߼�ֵ-1)��
% tempv = zeros(size(MT));
% for i = 1:length(MT)
%     tempv(i) = MT(i)./max(MT(1:i))-1;
% end
% value.zuidahuiche = min(tempv);
% %ƽ���Ƿ� ���㹫ʽΪ���˻��������ƽ��ֵ��
% value.pingjunzhangfu = mean(MT(2:end)./MT(1:end-1)/Wnum-1)/Wnum;
% %���Ǹ���-���㹫ʽΪ���������� / �ز⽻����������
% value.shangzhanggailv = length(find(diff(MT)>0))/(length(MT)-1);
% %��������������� ���������ǲ�������ӯ�������������
% testv1 = find(diff(MT)>0);
% %�ҳ��������ǵ������
% testv2 = diff(testv1);
% %��������Ų�Ϊ1�����ߴ���1
% testv3 = [(1:size(testv2))',testv2];
% testv3(eq(testv2,1),:)=[];
% value.maxv1 = max(diff(testv3(:,1)));
% %��������µ�����
% testv1 = find(diff(MT)<0);
% %�ҳ��µ��������
% testv2 = diff(testv1);
% %��������Ų�Ϊ1�����ߴ���1
% testv3 = [(1:size(testv2))',testv2];
% testv3(eq(testv2,1),:)=[];
% value.maxv2 = max(diff(testv3(:,1)));
% %���������������Ƿ����������ʲ�һ���ڵ���ñ��֡�
% value.maxdanqizhangfu = max(MT(2:end)./MT(1:end-1)-1);
% %������ڵ���
% value.maxdanqidiefu = min(MT(2:end)./MT(1:end-1)-1);
% %���沨����
% value.shouyibodonglv = std(MT(2:end)./MT(1:end-1)-1)*sqrt(52/Wnum);%�껯��׼��
