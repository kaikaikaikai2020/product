%{
Ҫ��
1 ƽ����ִ��
2 asia us ��ֵ����0.5������ִ�и��Ե�

3 ÿ�ζ�������ƽ��

�޸�1 ƽ��ǰ���ж��Ƿ����ƽ���Ƿ񴥷���
�޸�1 ƽ���ź�ΪԽ������
�ź�2-ģʽ2
%}

clear
close all
load data20201027.mat
symbol_pair = {'INFO','INFY','IRN_1M';'TT2330','TSM','NTN_1m'};
fee0 = [];
fee0.IN=[1/10000,1/10000];
fee0.TW = [3/10000,33/10000];
fee0.US=[1/10000,1/10000];
fees = {'IN','US';'TW','US'};

T1 = size(symbol_pair,1);
W0=20;
cut_v=1.5;
W = [W0-1,0];
f_str1 = '%s:����:���%0.4f��������ֵ%0.4f��������ֵ%0.4f %s';
f_str2 = '%s:����ƽ�����м�~�ǣ�%0.2f,US��%0.2f��ƽ�ּ�~���ޣ�%0.2f,US��%0.2f�����棺��%0.4f��US%0.4f';
f_str3 = '%s:����ƽ�����м�~�ǣ�%0.2f��ƽ�ּ�~���ޣ�%0.2f,���棺��%0.4f';
for i0 =1:1%T1
    
    recorder = strAdd();
    %data loading
    sym1 = symbol_pair{i0,1};
    sym2 = symbol_pair{i0,2};
    sym3 = symbol_pair{i0,3};
    title_str = sprintf('%s-%s',sym1,sym2);
    fn = sprintf('ADRƽ���ź�2ģʽ2-ִ�вֵ�%s_��ֵ%0.2f.xlsx',title_str,cut_v);
    
    x1 = data.(sym1);
    x2 = data.(sym2);
    x3 = data.(sym3);
    
    sub_fee1 = fee0.(fees{i0,1});
    sub_fee2 = fee0.(fees{i0,2});    
    sub_fee_com = [sub_fee1',sub_fee2']; 
    
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
    state = zeros(T2,1);%��¼״̬ 0 ���� 1 ���޶� 2 ���޿�
    for i = W0:T2
        %����ͳ��
        sub_pos1 = pos(i-1,:);%���ղ�λ�۸�
        sub_pos2 = [x1.data(i,2),x2.data(i,2)];%���տ��̼۸�
        sub_pos3 = x1.data(i,1);%�����������̼۸�
        sub_state = state(i-1);
        %�ж��Ƿ���Ҫƽ
        %����ƽ������������ƽ�� 
        q1 = eq(state(i-1),1) && boll(i-1)>=H(i-1) && boll(i)<=H(i);
        q2 = eq(state(i-1),-1) && boll(i-1)<=H(i-1) && boll(i)>=H(i);
        if q1 || q2
            if any(~eq(sub_pos1,0))
                r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                %��¼����
                y(i,:) = y(i,:)+r;
                %ƽ�֣���¼״̬
                sub_state = 0;
                sub_pos1(:) = 0;
            end 
            mark=1;
        else
            mark=0;
        end
        
        if boll(i)>H(i)
            if eq(sub_state,1)
                recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'�������ޣ��ֲֲ�ִ��'));
                %����
                pos(i,:) = sub_pos1;
                state(i) = sub_state;
            else
                recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
                %���Ͽ�������������
                pos(i,1)=sub_pos2(1);
                %����ͳ��
                temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1);
                if temp_r>0 %�����������0
                    y(i,1) = y(i,1)+temp_r; %��¼����
                    recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                    pos(i,1) = 0; %ƽ��
                    state(i) = 0; %״̬����
                else
                    pos(i,2) = -sub_pos2(2); %����������Ʊ
                    state(i) = 1;
                end
            end
        elseif boll(i)<L(i)
            if eq(sub_state,-1)
                 recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'�������ޣ��ֲֲ�ִ��'));
                %����
                pos(i,:) = sub_pos1;
                state(i) = sub_state;
            else
                recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
                %���Ͽ�������������
                pos(i,1)=-sub_pos2(1);
                %����ͳ��
                temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1); % �����������0
                if temp_r>0
                    y(i,1) = y(i,1)+temp_r;
                    recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                    pos(i,1) = 0;
                    state(i) = 0; %״̬����
                else
                    pos(i,2) = sub_pos2(2);
                    state(i) = -1; %״̬����
                end
            end
        else
            %����
            pos(i,:) = sub_pos1;
            state(i) = sub_state;
        end
    end
    y = y./2;
    %y1 = cumprod(1+sum(y,2));
    yc = [cumprod(1+y),cumprod(1+sum(y,2))];
    
    h = figure_S53(yc,tref,title_str,1);
    legend({sym1,sym2,[sym1,'-',sym2]},'NumColumns',3,'location','best');
    xlswrite(fn,recorder.str1);
end


