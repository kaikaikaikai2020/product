%{
1、信号1模式1   平仓信号为穿过bolling的上限或下限，如果位于bolling带的两端，如果有仓位不论仓位的方向都平仓再建仓，如果没有仓位直接建仓；
2、信号1模式2    平仓信号为穿过bolning的上限或下限，如果位于bolling带的两端，如果有仓位且方向相同，继续持有；
3、信号2模式1   平仓信号为穿过bolling的中线，如果位于bolling带的两端，如果有仓位不论仓位的方向都平仓再建仓，如果没有仓位直接建仓；
4、信号2模式2    平仓信号为穿过bolning的中线，如果位于bolling带的两端，如果有仓位且方向相同，继续持有；
%}
clear
close all
%参数设置
fee0 = [];
fee0.IN=[1/10000,1/10000];
fee0.TT = [3/10000,33/10000];
fee0.US=[1/10000,1/10000];
%载入数据
data = updatedataZX02();
%计算窗口
W0=20;
%保留结果点数
num_r = 60;
%阈值
for cut_v =[1.5, 2]
    YC = [];
    YC_info = [];
    H=cell(4,1);
    %rec = cell(4,1);
    rec = [];
    for method_num = 1:4

        if eq(method_num,1)
            method_id = 'signal1-model1';
        elseif eq(method_num,2)
            method_id = 'signal1-model2';
        elseif eq(method_num,3)
            method_id = 'signal2-model1';
        else
            method_id = 'signal2-model2';
        end

        T = length(data);
        re = cell(T,1);
        leg_str = re;
        re2= [];
        for i = 1:T    
            tempdata = data{i}; 
            sub_cp = tempdata.c_p;
            sym1 = tempdata.index{1};
            sym2 = tempdata.index{2};
            sym3 = tempdata.index{3};    
            
            x1 = tempdata.x1;
            x2 = tempdata.x2;
            x3 = tempdata.x3;    
            sub_fee1 = fee0.(tempdata.p_v{1});
            sub_fee2 = fee0.(tempdata.p_v{2});
            if eq(method_num,1)
                [tref,yc,y,recorder,info] =ADR_method.sig1mod1(x1,x2,x3,sub_fee1,sub_fee2,W0,cut_v,sub_cp);    
            elseif eq(method_num,2)
                [tref,yc,y,recorder,info] =ADR_method.sig1mod2(x1,x2,x3,sub_fee1,sub_fee2,W0,cut_v,sub_cp);    
            elseif eq(method_num,3)
                [tref,yc,y,recorder,info] =ADR_method.sig2mod1(x1,x2,x3,sub_fee1,sub_fee2,W0,cut_v,sub_cp);   
            else
                [tref,yc,y,recorder,info] =ADR_method.sig2mod2(x1,x2,x3,sub_fee1,sub_fee2,W0,cut_v,sub_cp);
            end

            re{i} = yc(:,end);
            leg_str{i} = sprintf('%s-%s',sym1,sym2);
            
            temp1 = info(end-num_r:end,:);
            temp2 = cell(size(temp1(:,1:2)));
            temp2(:,1) = {method_id};
            temp2(:,2) = {sprintf('%s-%s',sym1,sym2)};
            temp2 = cell2table(temp2,'VariableNames',{'methodID','pairName'});
            %rec{method_num} = [temp2,temp1];
            %rec = cat(1,rec,[temp2,temp1]);
            re2 =cat(1,re2, [temp2,temp1]);
        end

        yc=[re{:}];
        h1 = figure_S53(yc,tref,method_id,1);
        legend(leg_str,'Location','best');

        temp1 = yc(end-num_r:end,:);
        temp1 = bsxfun(@rdivide,temp1,temp1(1,:));
        h2 = figure_S53(temp1,tref(end-num_r:end),method_id,1);
        legend(leg_str,'Location','best');

        YC = cat(1,YC,re);
        YC_info = cat(1,YC_info,cellfun(@(x) sprintf('%s%s-',method_id,x), leg_str,'UniformOutput',false));
        H{method_num} = [h1,h2];

%         temp1 = info(end-num_r:end,:);
%         temp2 = cell(size(temp1(:,1)));
%         temp2(:) = {method_id};
%         temp2 = cell2table(temp2,'VariableNames',{'methodID'});
%         %rec{method_num} = [temp2,temp1];
%         rec = cat(1,rec,[temp2,temp1]);
          rec = cat(1,rec,re2);
    end
    H = [H{:}];
    sta_re = curve_static_batch(YC,YC_info);
    pn = fullfile(pwd,'计算结果');
    if ~exist(pn,'dir')
        mkdir(pn)
    end
    key_str = sprintf('支线任务2-ADR%s',datestr(now,'yyyymmdd'));
    fn1= fullfile(pn,sprintf('%s-%0.2f.doc',key_str,cut_v));
    obj = wordcom(fn1);
    for i = 1:length(H)
        obj.pasteFigure(H(i));
    end
    CloseWord(obj)

    fn2 = fullfile(pn,sprintf('%s曲线统计-%0.2f.csv',key_str,cut_v));
    writetable(cell2table(sta_re),fn2);

    fn3 = fullfile(pn,sprintf('%s操作明细-%0.2f.csv',key_str,cut_v));
    writetable(rec,fn3);
end

command = 'git add --all';
[status,cmdout] = dos(command)
command = 'git commit -a -m xxx';
[status,cmdout] = dos(command)
command = 'git push https://github.com/kaikaikaikai2020/test.git';
[status,cmdout] = dos(command)

pause(30);
command = 'git add --all';
[status,cmdout] = dos(command)
command = 'git commit -a -m xxx';
[status,cmdout] = dos(command)
command = 'git push https://github.com/kaikaikaikai2020/test.git';
[status,cmdout] = dos(command)

