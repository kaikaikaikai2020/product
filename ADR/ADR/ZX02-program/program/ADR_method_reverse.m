%update 20201231 Asia 和 US的收盘价格写到文件里面
%ntn+1m和irn+1m
%在原来ADR的基础上亚洲股和美股掉个，注意补充参数x3的对应关系

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
            %亚洲买入,美国买入
            y  = zeros(T2,2);
            f_str1 = '%s:开盘:溢价%0.4f，上限阈值%0.4f，下限阈值%0.4f %s';
            f_str2 = '%s:开盘平，持有价~America：%0.2f,亚洲：%0.2f，平仓价~America：%0.2f,亚洲：%0.2f，收益：America%0.4f，亚洲%0.4f';
            f_str3 = '%s:收盘平，持有价~America：%0.2f，平仓价~America：%0.2f,收益：America%0.4f';
            recorder = strAdd();
            oper = cell(T2,1);
            for i = W0:T2
                %开盘统计
                sub_pos1 = pos(i-1,:);%昨日仓位价格
                sub_pos2 = [x1.data(i,2),x2.data(i,2)];%今日开盘价格
                sub_pos3 = x1.data(i,1);%今日亚洲收盘价格                
                sub_recorder = strAdd();
                if boll(i)>H(i)
                    recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'超过上限'));
                    sub_recorder.A('超上限');
                    %符合平仓条件，开盘平仓            
                    if any(~eq(sub_pos1,0))
                        sub_recorder.A('开盘平');
                        r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                        recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                         %记录收益
                        y(i,:) = y(i,:)+r;
                        %sub_pos1(:) = 0;
                    end
                    %符合开仓条件，开仓
                    pos(i,1)=sub_pos2(1);
                    %收盘统计
                    temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1);
                    if temp_r>0 %收盘收益大于0
                        sub_recorder.A('收盘平');
                        y(i,1) = y(i,1)+temp_r; %记录收益
                        recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                        pos(i,1) = 0; %平仓
                    else
                        sub_recorder.A('收盘加仓亚股');
                        pos(i,2) = -sub_pos2(2); %买入美国股票
                    end
                elseif boll(i)<L(i)
                    recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'跌过下限'));
                    sub_recorder.A('跌下限');
                    %符合平仓条件，开盘平仓 
                    if any(~eq(sub_pos1,0))
                        sub_recorder.A('开盘平');
                        r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                        recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                        y(i,:)=y(i,:)+r;
                        %sub_pos1(:) = 0;
                    end
                    %符合开仓条件，开仓
                    pos(i,1)=-sub_pos2(1);
                    %收盘统计
                    temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1); % 收盘收益大于0
                    if temp_r>0
                        sub_recorder.A('收盘平');
                        y(i,1) = y(i,1)+temp_r;
                        recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                        pos(i,1) = 0;
                    else
                        sub_recorder.A('收盘加仓亚股');
                        pos(i,2) = sub_pos2(2);
                    end
                else
                    %保持
                    pos(i,:) = sub_pos1;
                end
                if ~isempty(sub_recorder.str1)>0
                    oper{i} = strjoin(sub_recorder.str1,'\');
                end
            end
            y = y./2;
            %y1 = cumprod(1+sum(y,2));
            yc = [cumprod(1+y),cumprod(1+sum(y,2))];
            %具体信息
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
            %亚洲买入,美国买入
            y  = zeros(T2,2);
            f_str1 = '%s:开盘:溢价%0.4f，上限阈值%0.4f，下限阈值%0.4f %s';
            f_str2 = '%s:开盘平，持有价~America：%0.2f,亚洲：%0.2f，平仓价~America：%0.2f,亚洲：%0.2f，收益：America%0.4f，亚洲%0.4f';
            f_str3 = '%s:收盘平，持有价~America：%0.2f，平仓价~America：%0.2f,收益：America%0.4f';
            
            recorder = strAdd();
            state = zeros(T2,1);%记录状态 0 正常 1 亚洲多 2 亚洲空
            oper = cell(T2,1);
            for i = W0:T2
                %开盘统计
                sub_pos1 = pos(i-1,:);%昨日仓位价格
                sub_pos2 = [x1.data(i,2),x2.data(i,2)];%今日开盘价格
                sub_pos3 = x1.data(i,1);%今日亚洲收盘价格
                sub_recorder = strAdd();
                if boll(i)>H(i)
                    if eq(state(i-1),1)
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'超过上限，持仓不执行'));
                        sub_recorder.A('超上限 持仓不执行');
                        %保持
                        pos(i,:) = pos(i-1,:);
                        state(i) = state(i-1);
                    else
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'超过上限'));
                        sub_recorder.A('超上限');
                        %符合平仓条件，开盘平仓                
                        if any(~eq(sub_pos1,0))
                            sub_recorder.A('开盘平');
                            r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                            recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                            %记录收益
                            y(i,:) = y(i,:)+r;
                            %平仓，记录状态
                            state(i) = 0;
                        end                
                        %符合开仓条件，开仓
                        pos(i,1)=sub_pos2(1);
                        %收盘统计
                        temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1);
                        if temp_r>0 %收盘收益大于0
                            sub_recorder.A('收盘平');
                            y(i,1) = y(i,1)+temp_r; %记录收益
                            recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                            pos(i,1) = 0; %平仓
                            state(i) = 0; %状态归零
                        else
                            sub_recorder.A('收盘加仓亚股');
                            pos(i,2) = -sub_pos2(2); %买入美国股票
                            state(i) = 1;
                        end
                    end
                elseif boll(i)<L(i)
                    if eq(state(i-1),-1)
                         recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'跌过下限，持仓不执行'));
                         sub_recorder.A('跌下限 持仓未执行');
                        %保持
                        pos(i,:) = pos(i-1,:);
                        state(i) = state(i-1);
                    else
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'跌过下限'));
                        sub_recorder.A('跌下限');
                        %符合平仓条件，开盘平仓
                        if any(~eq(sub_pos1,0))
                            sub_recorder.A('开盘平');
                            r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                            y(i,:)=y(i,:)+r;
                            recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                        end
                        %符合开仓条件，开仓
                        pos(i,1)=-sub_pos2(1);
                        %收盘统计
                        temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1); % 收盘收益大于0
                        if temp_r>0
                            sub_recorder.A('收盘平');
                            y(i,1) = y(i,1)+temp_r;
                            recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                            pos(i,1) = 0;
                            state(i) = 0; %状态归零
                        else
                            sub_recorder.A('收盘加仓亚股');
                            pos(i,2) = sub_pos2(2);
                            state(i) = -1; %状态归零
                        end
                    end
                else
                    %保持
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
            %亚洲买入,美国买入
            y  = zeros(T2,2);
            f_str1 = '%s:开盘:溢价%0.4f，上限阈值%0.4f，下限阈值%0.4f %s';
            f_str2 = '%s:开盘平，持有价~America：%0.2f,亚洲：%0.2f，平仓价~America：%0.2f,亚洲：%0.2f，收益：America%0.4f，亚洲%0.4f';
            f_str3 = '%s:收盘平，持有价~America：%0.2f，平仓价~America：%0.2f,收益：America%0.4f';
            
            recorder = strAdd();
            state = zeros(T2,1);%记录状态 0 正常 1 亚洲多 2 亚洲空
            oper = cell(T2,1);
            for i = W0:T2
                %开盘统计
                sub_pos1 = pos(i-1,:);%昨日仓位价格
                sub_pos2 = [x1.data(i,2),x2.data(i,2)];%今日开盘价格
                sub_pos3 = x1.data(i,1);%今日亚洲收盘价格
                sub_state = state(i-1);
                sub_recorder = strAdd();

                %判断是否需要平
                %符合平仓条件，开盘平仓 
                q1 = eq(sub_state,1) && boll(i)<H(i);
                q2 = eq(sub_state,-1) &&  boll(i)>H(i);
                if q1 || q2
                    if q1
                        sub_recorder.A('下穿中线');
                    else
                        sub_recorder.A('上穿中线');
                    end
                    if any(~eq(sub_pos1,0))
                        sub_recorder.A('开盘平');
                        r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                        recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                        %记录收益
                        y(i,:) = y(i,:)+r;
                        %平仓，记录状态
                        sub_state = 0;
                        sub_pos1(:) = 0;
                    end 
                    mark=1;
                else
                    mark=0;
                end

                if boll(i)>H(i)
                    recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'超过上限'));
                    sub_recorder.A('超上限');
                    %是否有仓位
                    if any(~eq(sub_pos1,0)) && eq(mark,0)
                        sub_recorder.A('开盘平');
                        r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                        recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                        %记录收益
                        y(i,:) = y(i,:)+r;
                        %平仓，记录状态
                        state(i) = 0;
                    end  
                    %符合开仓条件，开仓
                    pos(i,1)=sub_pos2(1);
                    %收盘统计
                    temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1);
                    if temp_r>0 %收盘收益大于0
                        sub_recorder.A('收盘平');
                        y(i,1) = y(i,1)+temp_r; %记录收益
                        recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                        pos(i,1) = 0; %平仓
                        state(i) = 0; %状态归零
                    else
                        sub_recorder.A('收盘加仓亚股');
                        pos(i,2) = -sub_pos2(2); %买入美国股票
                        state(i) = 1;
                    end
                elseif boll(i)<L(i)
                    recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'跌过下限'));
                    sub_recorder.A('跌下限');
                    %判断是否有仓位
                    if any(~eq(sub_pos1,0)) && eq(mark,0)
                        sub_recorder.A('开盘平');
                        r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                        y(i,:)=y(i,:)+r;
                        recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                        state(i) = 0;
                    end
                    %符合开仓条件，开仓
                    pos(i,1)=-sub_pos2(1);
                    %收盘统计
                    temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1); % 收盘收益大于0
                    if temp_r>0
                        sub_recorder.A('收盘平');
                        y(i,1) = y(i,1)+temp_r;
                        recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                        pos(i,1) = 0;
                        state(i) = 0; %状态归零
                    else
                        sub_recorder.A('收盘加仓亚股');
                        pos(i,2) = sub_pos2(2);
                        state(i) = -1; %状态归零
                    end
                else
                    %保持
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
            %具体信息
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
            %亚洲买入,美国买入
            y  = zeros(T2,2);
            f_str1 = '%s:开盘:溢价%0.4f，上限阈值%0.4f，下限阈值%0.4f %s';
            f_str2 = '%s:开盘平，持有价~America：%0.2f,亚洲：%0.2f，平仓价~America：%0.2f,亚洲：%0.2f，收益：America%0.4f，亚洲%0.4f';
            f_str3 = '%s:收盘平，持有价~America：%0.2f，平仓价~America：%0.2f,收益：America%0.4f';

            recorder = strAdd();
            state = zeros(T2,1);%记录状态 0 正常 1 亚洲多 2 亚洲空
            oper = cell(T2,1);
            for i = W0:T2
                %开盘统计
                sub_pos1 = pos(i-1,:);%昨日仓位价格
                sub_pos2 = [x1.data(i,2),x2.data(i,2)];%今日开盘价格
                sub_pos3 = x1.data(i,1);%今日亚洲收盘价格
                sub_state = state(i-1);
                sub_recorder = strAdd();
                %判断是否需要平
                %符合平仓条件，开盘平仓 
                q1 = eq(sub_state,1) && boll(i)<H(i);
                q2 = eq(sub_state,-1) &&  boll(i)>H(i);
                if q1 || q2
                    if q1
                        sub_recorder.A('下穿中线');
                    else
                        sub_recorder.A('上穿中线');
                    end
                    if any(~eq(sub_pos1,0))
                        sub_recorder.A('开盘平');
                        r = get_ret(sub_pos1,sub_pos2,sub_fee_com);
                        recorder.A(sprintf(f_str2,tref{i},sub_pos1,sub_pos2,r));
                        %记录收益
                        y(i,:) = y(i,:)+r;
                        %平仓，记录状态
                        sub_state = 0;
                        sub_pos1(:) = 0;
                    end 
                end

                if boll(i)>H(i)
                    if eq(sub_state,1)
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'超过上限，持仓不执行'));
                        sub_recorder.A('超上限 持仓不执行');
                        %保持
                        pos(i,:) = sub_pos1;
                        state(i) = sub_state;
                    else
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'超过上限'));
                        sub_recorder.A('超上限');
                        %符合开仓条件，开仓
                        pos(i,1)=sub_pos2(1);
                        %收盘统计
                        temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1);
                        if temp_r>0 %收盘收益大于0
                            sub_recorder.A('收盘平');
                            y(i,1) = y(i,1)+temp_r; %记录收益
                            recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                            pos(i,1) = 0; %平仓
                            state(i) = 0; %状态归零
                        else
                            sub_recorder.A('收盘加仓亚股');
                            pos(i,2) = -sub_pos2(2); %买入美国股票
                            state(i) = 1;
                        end
                    end
                elseif boll(i)<L(i)
                    if eq(sub_state,-1)
                         recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'跌过下限，持仓不执行'));
                         sub_recorder.A('跌下限 持仓不执行');
                        %保持
                        pos(i,:) = sub_pos1;
                        state(i) = sub_state;
                    else
                        recorder.A(sprintf(f_str1,tref{i},boll(i),H(i),L(i),'跌过下限'));
                        sub_recorder.A('跌下限');
                        %符合开仓条件，开仓
                        pos(i,1)=-sub_pos2(1);
                        %收盘统计
                        temp_r = get_ret0(pos(i,1),sub_pos3,sub_fee1); % 收盘收益大于0
                        if temp_r>0
                            sub_recorder.A('收盘平');
                            y(i,1) = y(i,1)+temp_r;
                            recorder.A(sprintf(f_str3,tref{i},pos(i,1),sub_pos3,temp_r));
                            pos(i,1) = 0;
                            state(i) = 0; %状态归零
                        else
                            sub_recorder.A('收盘加仓亚股');
                            pos(i,2) = sub_pos2(2);
                            state(i) = -1; %状态归零
                        end
                    end
                else
                    %保持
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
            %具体信息
%             info = [tref,num2cell([boll,H,L,M,pos]),oper];
%             info = cell2table(info,'VariableNames',{'date','premium','HighBond','LowBond','meanLine','asiaPos','usPos','operation'});
            info = [tref,num2cell([boll,H,L,M,HP,LP,pos,x1.data(:,1:2),x2.data(:,1:2),x3.data(:,1:2)]),oper];
            info = cell2table(info,'VariableNames',{'date','premium','HighBond','LowBond','meanLine','Hbond_price','Lbond_price','USPos','asiaPos',...
                'USOP','USCL','asiaOP','asiaCL','r_exOP','r_exCL','operation'});
            
        end
    end
end