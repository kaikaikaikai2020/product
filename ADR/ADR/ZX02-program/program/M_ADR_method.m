%{
Ҫ��
1 ƽ����ִ��
2 asia us ��ֵ����0.5������ִ�и��Ե�
%}

clear
close all
load data20201027.mat
symbol_pair = {'INFO','INFY','IRN_1M';'TT2330','TSM','NTN_1m'};
T1 = size(symbol_pair,1);
W0=20;
cut_v=1.5;
W = [W0-1,0];
f_str1 = '%s:����:���%0.4f��������ֵ%0.4f��������ֵ%0.4f %s';
f_str2 = '%s:����ƽ�����м�~�ǣ�%0.2f,US��%0.2f��ƽ�ּ�~���ޣ�%0.2f,US��%0.2f�����棺��%0.4f��US%0.4f';
f_str3 = '%s:����ƽ�����м�~�ǣ�%0.2f��ƽ�ּ�~���ޣ�%0.2f,���棺��%0.4f';
for i0 =1:1%T1
    obj = strAdd();
    %data loading
    sym1 = symbol_pair{i0,1};
    sym2 = symbol_pair{i0,2};
    sym3 = symbol_pair{i0,3};
    title_str = sprintf('%s-%s',sym1,sym2);
    x1 = data.(sym1);
    x2 = data.(sym2);
    x3 = data.(sym3);
    tref = x1.tref;
    open_asia = x1.data(:,2);
    close_us = x2.data(:,1);
    p = x3.data(:,2);
    %signal
    boll = zeros(size(close_us));
    boll(2:end) = close_us(1:end-1).*p(2:end)./open_asia(2:end)-1;
    M = movmean(boll,W);
    S = movstd(boll,W);
    L = M-S*cut_v;
    H = M+S*cut_v;

    T2 =length(boll);
    pos=zeros(T2,2);
    %��������,��������
    y  = zeros(T2,2);
    for i = W0:T2
        %����ͳ��
        sub_pos1 = pos(i-1,:);%���ղ�λ�۸�
        sub_pos2 = [x1.data(i,2),x2.data(i,2)];%���տ��̼۸�
        sub_pos3 = x1.data(i,1);%�����������̼۸�
        
        if boll(i)>H(i)
            obj.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
            %����ƽ������������ƽ��
            r = get_ret(sub_pos1,sub_pos2);
            if any(~eq(sub_pos1,0))
                obj.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
            end
            %��¼����
            y(i,:) = y(i,:)+r;
            %���Ͽ�������������
            pos(i,1)=sub_pos2(1);
            %����ͳ��
            temp_r = get_ret0(pos(i,1),sub_pos3);
            if temp_r>0 %�����������0
                y(i,1) = y(i,1)+temp_r; %��¼����
                obj.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                pos(i,1) = 0; %ƽ��
            else
                pos(i,2) = -sub_pos2(2); %����������Ʊ
            end
        elseif boll(i)<L(i)
            obj.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
            %����ƽ������������ƽ��     
            r = get_ret(sub_pos1,sub_pos2);
            if any(~eq(sub_pos1,0))
                obj.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
            end
            y(i,:)=y(i,:)+r;
            %���Ͽ�������������
            pos(i,1)=-sub_pos2(1);
            %����ͳ��
            temp_r = get_ret0(pos(i,1),sub_pos3); % �����������0
            if temp_r>0
                y(i,1) = y(i,1)+temp_r;
                obj.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                pos(i,1) = 0;
            else
                pos(i,2) = sub_pos2(2);
            end
        else
            %����
            pos(i,:) = pos(i-1,:);
        end
    end
    %y1 = cumprod(1+sum(y,2));
    y = [0.5,0.5].*cumprod(1+y);
    
    h = figure_S53([y,sum(y,2)],tref,title_str,1);
    legend({sym1,sym2,[sym1,'-',sym2]},'NumColumns',3,'location','best');
    xlswrite(sprintf('ִ�вֵ�%s_��ֵ%0.2f.xlsx',title_str,cut_v),obj.str1);
end


