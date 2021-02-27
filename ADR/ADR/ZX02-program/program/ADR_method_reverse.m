%update 20201231 Asia �� US�����̼۸�д���ļ�����
%ntn+1m��irn+1m
%��ԭ��ADR�Ļ��������޹ɺ����ɵ�����ע�ⲹ�����x3�Ķ�Ӧ��ϵ

classdef ADR_method_reverse
    
    methods(Static)        
        function  [tref,yc,y,recorder,info] =sig1mod1(x1,x2,x3,sub_fee1,sub_fee2,W0,cut_v,sub_cp)
            if nargin < 4
                sub_fee1 = [1,1]/10000;
            end
            if nargin < 5
                sub_fee2 = sub_fee1;
            end
            if nargin < 6
                W0 = 20;
            end
            if nargin < 7
                cut_v=1.5;
            end
            if nargin <8
                sub_cp=1;
            end
            W = [W0-1,0];
            sub_fee_com = [sub_fee1',sub_fee2'];  
            
            tref = x1.tref;
            open_asia = x1.data(:,2);
            close_us = x2.data(:,1);
            p = x3.data(:,2);
            %signal
            boll = zeros(size(close_us));
            boll(2:end) = close_us(1:end-1).*p(2:end)./open_asia(2:end)/sub_cp-1;
            M = movmean(boll,W);
            S = movstd(boll,W);
            L = M-S*cut_v;
            H = M+S*cut_v;
            
            HP=zeros(size(H));
            LP=zeros(size(L));
            HP(2:end) = close_us(1:end-1).*p(2:end)./sub_cp./(H(2:end)+1);
            LP(2:end) = close_us(1:end-1).*p(2:end)./sub_cp./(L(2:end)+1);
            T2 =length(boll);
            pos=zeros(T2,2);
            %��������,��������
            y  = zeros(T2,2);
            f_str1 = '%s:����:���%0.4f��������ֵ%0.4f��������ֵ%0.4f %s';
            f_str2 = '%s:����ƽ�����м�~America��%0.2f,���ޣ�%0.2f��ƽ�ּ�~America��%0.2f,���ޣ�%0.2f�����棺America%0.4f������%0.4f';
            f_str3 = '%s:����ƽ�����м�~America��%0.2f��ƽ�ּ�~America��%0.2f,���棺America%0.4f';
            recorder = strAdd();
            oper = cell(T2,1);
            for i = W0:T2
                %����ͳ��
                sub_pos1 = pos(i-1,:);%���ղ�λ�۸�
                sub_pos2 = [x1.data(i,2),x2.data(i,2)];%���տ��̼۸�
                sub_pos3 = x1.data(i,1);%�����������̼۸�                
                sub_recorder = strAdd();
                if boll(i)>H(i)
                    recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
                    sub_recorder.A('������');
                    %����ƽ������������ƽ��            
                    if any(~eq(sub_pos1,0))
                        sub_recorder.A('����ƽ');
                        r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                        recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                         %��¼����
                        y(i,:) = y(i,:)+r;
                        %sub_pos1(:) = 0;
                    end
                    %���Ͽ�������������
                    pos(i,1)=sub_pos2(1);
                    %����ͳ��
                    temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1);
                    if temp_r>0 %�����������0
                        sub_recorder.A('����ƽ');
                        y(i,1) = y(i,1)+temp_r; %��¼����
                        recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                        pos(i,1) = 0; %ƽ��
                    else
                        sub_recorder.A('���̼Ӳ��ǹ�');
                        pos(i,2) = -sub_pos2(2); %����������Ʊ
                    end
                elseif boll(i)<L(i)
                    recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
                    sub_recorder.A('������');
                    %����ƽ������������ƽ�� 
                    if any(~eq(sub_pos1,0))
                        sub_recorder.A('����ƽ');
                        r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                        recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                        y(i,:)=y(i,:)+r;
                        %sub_pos1(:) = 0;
                    end
                    %���Ͽ�������������
                    pos(i,1)=-sub_pos2(1);
                    %����ͳ��
                    temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1); % �����������0
                    if temp_r>0
                        sub_recorder.A('����ƽ');
                        y(i,1) = y(i,1)+temp_r;
                        recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                        pos(i,1) = 0;
                    else
                        sub_recorder.A('���̼Ӳ��ǹ�');
                        pos(i,2) = sub_pos2(2);
                    end
                else
                    %����
                    pos(i,:) = sub_pos1;
                end
                if ~isempty(sub_recorder.str1)>0
                    oper{i} = strjoin(sub_recorder.str1,'\');
                end
            end
            y = y./2;
            %y1 = cumprod(1+sum(y,2));
            yc = [cumprod(1+y),cumprod(1+sum(y,2))];
            %������Ϣ
            info = [tref,num2cell([boll,H,L,M,HP,LP,pos,x1.data(:,1:2),x2.data(:,1:2),x3.data(:,1:2)]),oper];
            info = cell2table(info,'VariableNames',{'date','premium','HighBond','LowBond','meanLine','Hbond_price','Lbond_price','USPos','asiaPos',...
                'USOP','USCL','asiaOP','asiaCL','r_exOP','r_exCL','operation'});
            
        end
        function  [tref,yc,y,recorder,info] =sig1mod2(x1,x2,x3,sub_fee1,sub_fee2,W0,cut_v,sub_cp)
            if nargin < 4
                sub_fee1 = [1,1]/10000;
            end
            if nargin < 5
                sub_fee2 = sub_fee1;
            end
            if nargin < 6
                W0 = 20;
            end
            if nargin < 7
                cut_v=1.5;
            end
            if nargin <8
                sub_cp=1;
            end
            W = [W0-1,0];
            sub_fee_com = [sub_fee1',sub_fee2'];  
            
            tref = x1.tref;
            open_asia = x1.data(:,2);
            close_us = x2.data(:,1);
            p = x3.data(:,2);
            %signal
            boll = zeros(size(close_us));
            boll(2:end) = close_us(1:end-1).*p(2:end)./open_asia(2:end)/sub_cp-1;
            M = movmean(boll,W);
            S = movstd(boll,W);
            L = M-S*cut_v;
            H = M+S*cut_v;
            
            HP=zeros(size(H));
            LP=zeros(size(L));
            HP(2:end) = close_us(1:end-1).*p(2:end)./sub_cp./(H(2:end)+1);
            LP(2:end) = close_us(1:end-1).*p(2:end)./sub_cp./(L(2:end)+1);
            
            T2 =length(boll);
            pos=zeros(T2,2);
            %��������,��������
            y  = zeros(T2,2);
            f_str1 = '%s:����:���%0.4f��������ֵ%0.4f��������ֵ%0.4f %s';
            f_str2 = '%s:����ƽ�����м�~America��%0.2f,���ޣ�%0.2f��ƽ�ּ�~America��%0.2f,���ޣ�%0.2f�����棺America%0.4f������%0.4f';
            f_str3 = '%s:����ƽ�����м�~America��%0.2f��ƽ�ּ�~America��%0.2f,���棺America%0.4f';
            
            recorder = strAdd();
            state = zeros(T2,1);%��¼״̬ 0 ���� 1 ���޶� 2 ���޿�
            oper = cell(T2,1);
            for i = W0:T2
                %����ͳ��
                sub_pos1 = pos(i-1,:);%���ղ�λ�۸�
                sub_pos2 = [x1.data(i,2),x2.data(i,2)];%���տ��̼۸�
                sub_pos3 = x1.data(i,1);%�����������̼۸�
                sub_recorder = strAdd();
                if boll(i)>H(i)
                    if eq(state(i-1),1)
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'�������ޣ��ֲֲ�ִ��'));
                        sub_recorder.A('������ �ֲֲ�ִ��');
                        %����
                        pos(i,:) = pos(i-1,:);
                        state(i) = state(i-1);
                    else
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
                        sub_recorder.A('������');
                        %����ƽ������������ƽ��                
                        if any(~eq(sub_pos1,0))
                            sub_recorder.A('����ƽ');
                            r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                            recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                            %��¼����
                            y(i,:) = y(i,:)+r;
                            %ƽ�֣���¼״̬
                            state(i) = 0;
                        end                
                        %���Ͽ�������������
                        pos(i,1)=sub_pos2(1);
                        %����ͳ��
                        temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1);
                        if temp_r>0 %�����������0
                            sub_recorder.A('����ƽ');
                            y(i,1) = y(i,1)+temp_r; %��¼����
                            recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                            pos(i,1) = 0; %ƽ��
                            state(i) = 0; %״̬����
                        else
                            sub_recorder.A('���̼Ӳ��ǹ�');
                            pos(i,2) = -sub_pos2(2); %����������Ʊ
                            state(i) = 1;
                        end
                    end
                elseif boll(i)<L(i)
                    if eq(state(i-1),-1)
                         recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'�������ޣ��ֲֲ�ִ��'));
                         sub_recorder.A('������ �ֲ�δִ��');
                        %����
                        pos(i,:) = pos(i-1,:);
                        state(i) = state(i-1);
                    else
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
                        sub_recorder.A('������');
                        %����ƽ������������ƽ��
                        if any(~eq(sub_pos1,0))
                            sub_recorder.A('����ƽ');
                            r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                            y(i,:)=y(i,:)+r;
                            recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                        end
                        %���Ͽ�������������
                        pos(i,1)=-sub_pos2(1);
                        %����ͳ��
                        temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1); % �����������0
                        if temp_r>0
                            sub_recorder.A('����ƽ');
                            y(i,1) = y(i,1)+temp_r;
                            recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                            pos(i,1) = 0;
                            state(i) = 0; %״̬����
                        else
                            sub_recorder.A('���̼Ӳ��ǹ�');
                            pos(i,2) = sub_pos2(2);
                            state(i) = -1; %״̬����
                        end
                    end
                else
                    %����
                    pos(i,:) = pos(i-1,:);
                    state(i) = state(i-1);
                end
                if ~isempty(sub_recorder.str1)>0
                    oper{i} = strjoin(sub_recorder.str1,'\');
                end
            end
            y = y./2;
            %y1 = cumprod(1+sum(y,2));
            yc = [cumprod(1+y),cumprod(1+sum(y,2))];
%             info = [tref,num2cell([boll,H,L,M,pos]),oper];
%             info = cell2table(info,'VariableNames',{'date','premium','HighBond','LowBond','meanLine','asiaPos','usPos','operation'});
            info = [tref,num2cell([boll,H,L,M,HP,LP,pos,x1.data(:,1:2),x2.data(:,1:2),x3.data(:,1:2)]),oper];
            info = cell2table(info,'VariableNames',{'date','premium','HighBond','LowBond','meanLine','Hbond_price','Lbond_price','USPos','asiaPos',...
                'USOP','USCL','asiaOP','asiaCL','r_exOP','r_exCL','operation'});
            
        end
        
        function  [tref,yc,y,recorder,info] =sig2mod1(x1,x2,x3,sub_fee1,sub_fee2,W0,cut_v,sub_cp)
            if nargin < 4
                sub_fee1 = [1,1]/10000;
            end
            if nargin < 5
                sub_fee2 = sub_fee1;
            end
            if nargin < 6
                W0 = 20;
            end
            if nargin < 7
                cut_v=1.5;
            end
            if nargin <8
                sub_cp=1;
            end
            W = [W0-1,0];
            sub_fee_com = [sub_fee1',sub_fee2'];  
            
            tref = x1.tref;
            open_asia = x1.data(:,2);
            close_us = x2.data(:,1);
            p = x3.data(:,2);
            %signal
            boll = zeros(size(close_us));
            boll(2:end) = close_us(1:end-1).*p(2:end)./open_asia(2:end)/sub_cp-1;
            M = movmean(boll,W);
            S = movstd(boll,W);
            L = M-S*cut_v;
            H = M+S*cut_v;
            
            HP=zeros(size(H));
            LP=zeros(size(L));
            HP(2:end) = close_us(1:end-1).*p(2:end)./sub_cp./(H(2:end)+1);
            LP(2:end) = close_us(1:end-1).*p(2:end)./sub_cp./(L(2:end)+1);
            
            T2 =length(boll);
            pos=zeros(T2,2);
            %��������,��������
            y  = zeros(T2,2);
            f_str1 = '%s:����:���%0.4f��������ֵ%0.4f��������ֵ%0.4f %s';
            f_str2 = '%s:����ƽ�����м�~America��%0.2f,���ޣ�%0.2f��ƽ�ּ�~America��%0.2f,���ޣ�%0.2f�����棺America%0.4f������%0.4f';
            f_str3 = '%s:����ƽ�����м�~America��%0.2f��ƽ�ּ�~America��%0.2f,���棺America%0.4f';
            
            recorder = strAdd();
            state = zeros(T2,1);%��¼״̬ 0 ���� 1 ���޶� 2 ���޿�
            oper = cell(T2,1);
            for i = W0:T2
                %����ͳ��
                sub_pos1 = pos(i-1,:);%���ղ�λ�۸�
                sub_pos2 = [x1.data(i,2),x2.data(i,2)];%���տ��̼۸�
                sub_pos3 = x1.data(i,1);%�����������̼۸�
                sub_state = state(i-1);
                sub_recorder = strAdd();

                %�ж��Ƿ���Ҫƽ
                %����ƽ������������ƽ�� 
                q1 = eq(sub_state,1) && boll(i)<H(i);
                q2 = eq(sub_state,-1) &&  boll(i)>H(i);
                if q1 || q2
                    if q1
                        sub_recorder.A('�´�����');
                    else
                        sub_recorder.A('�ϴ�����');
                    end
                    if any(~eq(sub_pos1,0))
                        sub_recorder.A('����ƽ');
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
                    recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
                    sub_recorder.A('������');
                    %�Ƿ��в�λ
                    if any(~eq(sub_pos1,0)) && eq(mark,0)
                        sub_recorder.A('����ƽ');
                        r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                        recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                        %��¼����
                        y(i,:) = y(i,:)+r;
                        %ƽ�֣���¼״̬
                        state(i) = 0;
                    end  
                    %���Ͽ�������������
                    pos(i,1)=sub_pos2(1);
                    %����ͳ��
                    temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1);
                    if temp_r>0 %�����������0
                        sub_recorder.A('����ƽ');
                        y(i,1) = y(i,1)+temp_r; %��¼����
                        recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                        pos(i,1) = 0; %ƽ��
                        state(i) = 0; %״̬����
                    else
                        sub_recorder.A('���̼Ӳ��ǹ�');
                        pos(i,2) = -sub_pos2(2); %����������Ʊ
                        state(i) = 1;
                    end
                elseif boll(i)<L(i)
                    recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
                    sub_recorder.A('������');
                    %�ж��Ƿ��в�λ
                    if any(~eq(sub_pos1,0)) && eq(mark,0)
                        sub_recorder.A('����ƽ');
                        r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                        y(i,:)=y(i,:)+r;
                        recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                        state(i) = 0;
                    end
                    %���Ͽ�������������
                    pos(i,1)=-sub_pos2(1);
                    %����ͳ��
                    temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1); % �����������0
                    if temp_r>0
                        sub_recorder.A('����ƽ');
                        y(i,1) = y(i,1)+temp_r;
                        recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                        pos(i,1) = 0;
                        state(i) = 0; %״̬����
                    else
                        sub_recorder.A('���̼Ӳ��ǹ�');
                        pos(i,2) = sub_pos2(2);
                        state(i) = -1; %״̬����
                    end
                else
                    %����
                    pos(i,:) = sub_pos1;
                    state(i) = sub_state;
                end
                if ~isempty(sub_recorder.str1)>0
                    oper{i} = strjoin(sub_recorder.str1,'\');
                end
            end
            y = y./2;
            %y1 = cumprod(1+sum(y,2));
            yc = [cumprod(1+y),cumprod(1+sum(y,2))];
            %������Ϣ
%             info = [tref,num2cell([boll,H,L,M,pos]),oper];
%             info = cell2table(info,'VariableNames',{'date','premium','HighBond','LowBond','meanLine','asiaPos','usPos','operation'});
            info = [tref,num2cell([boll,H,L,M,HP,LP,pos,x1.data(:,1:2),x2.data(:,1:2),x3.data(:,1:2)]),oper];
            info = cell2table(info,'VariableNames',{'date','premium','HighBond','LowBond','meanLine','Hbond_price','Lbond_price','USPos','asiaPos',...
                'USOP','USCL','asiaOP','asiaCL','r_exOP','r_exCL','operation'});
        end
        
        function  [tref,yc,y,recorder,info] =sig2mod2(x1,x2,x3,sub_fee1,sub_fee2,W0,cut_v,sub_cp)
            if nargin < 4
                sub_fee1 = [1,1]/10000;
            end
            if nargin < 5
                sub_fee2 = sub_fee1;
            end
            if nargin < 6
                W0 = 20;
            end
            if nargin < 7
                cut_v=1.5;
            end
            if nargin <8
                sub_cp=1;
            end
            W = [W0-1,0];
            sub_fee_com = [sub_fee1',sub_fee2'];  
            
            tref = x1.tref;
            open_asia = x1.data(:,2);
            close_us = x2.data(:,1);
            p = x3.data(:,2);
            %signal
            boll = zeros(size(close_us));
            boll(2:end) = close_us(1:end-1).*p(2:end)./open_asia(2:end)/sub_cp-1;
            M = movmean(boll,W);
            S = movstd(boll,W);
            L = M-S*cut_v;
            H = M+S*cut_v;
            
            HP=zeros(size(H));
            LP=zeros(size(L));
            HP(2:end) = close_us(1:end-1).*p(2:end)./sub_cp./(H(2:end)+1);
            LP(2:end) = close_us(1:end-1).*p(2:end)./sub_cp./(L(2:end)+1);
            T2 =length(boll);
            pos=zeros(T2,2);
            %��������,��������
            y  = zeros(T2,2);
            f_str1 = '%s:����:���%0.4f��������ֵ%0.4f��������ֵ%0.4f %s';
            f_str2 = '%s:����ƽ�����м�~America��%0.2f,���ޣ�%0.2f��ƽ�ּ�~America��%0.2f,���ޣ�%0.2f�����棺America%0.4f������%0.4f';
            f_str3 = '%s:����ƽ�����м�~America��%0.2f��ƽ�ּ�~America��%0.2f,���棺America%0.4f';

            recorder = strAdd();
            state = zeros(T2,1);%��¼״̬ 0 ���� 1 ���޶� 2 ���޿�
            oper = cell(T2,1);
            for i = W0:T2
                %����ͳ��
                sub_pos1 = pos(i-1,:);%���ղ�λ�۸�
                sub_pos2 = [x1.data(i,2),x2.data(i,2)];%���տ��̼۸�
                sub_pos3 = x1.data(i,1);%�����������̼۸�
                sub_state = state(i-1);
                sub_recorder = strAdd();
                %�ж��Ƿ���Ҫƽ
                %����ƽ������������ƽ�� 
                q1 = eq(sub_state,1) && boll(i)<H(i);
                q2 = eq(sub_state,-1) &&  boll(i)>H(i);
                if q1 || q2
                    if q1
                        sub_recorder.A('�´�����');
                    else
                        sub_recorder.A('�ϴ�����');
                    end
                    if any(~eq(sub_pos1,0))
                        sub_recorder.A('����ƽ');
                        r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                        recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                        %��¼����
                        y(i,:) = y(i,:)+r;
                        %ƽ�֣���¼״̬
                        sub_state = 0;
                        sub_pos1(:) = 0;
                    end 
                end

                if boll(i)>H(i)
                    if eq(sub_state,1)
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'�������ޣ��ֲֲ�ִ��'));
                        sub_recorder.A('������ �ֲֲ�ִ��');
                        %����
                        pos(i,:) = sub_pos1;
                        state(i) = sub_state;
                    else
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
                        sub_recorder.A('������');
                        %���Ͽ�������������
                        pos(i,1)=sub_pos2(1);
                        %����ͳ��
                        temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1);
                        if temp_r>0 %�����������0
                            sub_recorder.A('����ƽ');
                            y(i,1) = y(i,1)+temp_r; %��¼����
                            recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                            pos(i,1) = 0; %ƽ��
                            state(i) = 0; %״̬����
                        else
                            sub_recorder.A('���̼Ӳ��ǹ�');
                            pos(i,2) = -sub_pos2(2); %����������Ʊ
                            state(i) = 1;
                        end
                    end
                elseif boll(i)<L(i)
                    if eq(sub_state,-1)
                         recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'�������ޣ��ֲֲ�ִ��'));
                         sub_recorder.A('������ �ֲֲ�ִ��');
                        %����
                        pos(i,:) = sub_pos1;
                        state(i) = sub_state;
                    else
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'��������'));
                        sub_recorder.A('������');
                        %���Ͽ�������������
                        pos(i,1)=-sub_pos2(1);
                        %����ͳ��
                        temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1); % �����������0
                        if temp_r>0
                            sub_recorder.A('����ƽ');
                            y(i,1) = y(i,1)+temp_r;
                            recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                            pos(i,1) = 0;
                            state(i) = 0; %״̬����
                        else
                            sub_recorder.A('���̼Ӳ��ǹ�');
                            pos(i,2) = sub_pos2(2);
                            state(i) = -1; %״̬����
                        end
                    end
                else
                    %����
                    pos(i,:) = sub_pos1;
                    state(i) = sub_state;
                end
                if ~isempty(sub_recorder.str1)>0
                    oper{i} = strjoin(sub_recorder.str1,'\');
                end
            end
            y = y./2;
            %y1 = cumprod(1+sum(y,2));
            yc = [cumprod(1+y),cumprod(1+sum(y,2))];
            %������Ϣ
%             info = [tref,num2cell([boll,H,L,M,pos]),oper];
%             info = cell2table(info,'VariableNames',{'date','premium','HighBond','LowBond','meanLine','asiaPos','usPos','operation'});
            info = [tref,num2cell([boll,H,L,M,HP,LP,pos,x1.data(:,1:2),x2.data(:,1:2),x3.data(:,1:2)]),oper];
            info = cell2table(info,'VariableNames',{'date','premium','HighBond','LowBond','meanLine','Hbond_price','Lbond_price','USPos','asiaPos',...
                'USOP','USCL','asiaOP','asiaCL','r_exOP','r_exCL','operation'});
            
        end
    end
end