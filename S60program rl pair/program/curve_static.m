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
%1年化收益率
v_str{1} = '年化收益率';
v(1) = (y(end)/y(1))^(wid_year/length(y))-1;
sta_val.nh = v(1);
%2最大回撤
v_str{2} = '最大回撤';
v(2) = min(getdrawdown(y));
sta_val.drawdown = v(2);
%3-5持有一年最小\最大\平均平均收益
if length(y)>=wid_year
    v1 = y(wid_year:end)./y(1:end-wid_year+1)-1;
else
    v1 = nan;
end
v_str{3} = '持有1年最低收益';
v(3) = min(v1);
v_str{4} = '持有1年最高收益';
v(4) = max(v1);
v_str{5} = '持有1年平均收益';
v(5) = mean(v1);
sta_val.min_return = v(3);
sta_val.max_return = v(4);
sta_val.ave_return = v(5);
%6最大连续上涨数 衡量了我们策略连续盈利的最大数。
testv1 = find(diff(y)>0);
    %找出连续上涨的周序号
testv2 = diff(testv1);
    %连续的序号差为1，否者大于1
testv3 = [(1:size(testv2))',testv2];
testv3(eq(testv2,1),:)=[];
v_str{6} = '最大连续上涨天数';
if  size(testv3,1)>1
    v(6) = max(diff(testv3(:,1)));
else
    v(6) = 1;
end
sta_val.max_suc_up = v(6);
%7最大连续下跌数
testv1 = find(diff(y)<0);
    %找出下跌的序号
testv2 = diff(testv1);
    %连续的序号差为1，否者大于1
testv3 = [(1:size(testv2))',testv2];
testv3(eq(testv2,1),:)=[];
v_str{7} = '最大连续下跌天数';
if size(testv3,1)>1
    v(7) = max(diff(testv3(:,1)));
else
    v(7) = 1;
end
sta_val.max_suc_down = v(7);
%8收益波动率
    %std(每日收益率)*sqrt(252)
v_str{8} = '收益波动率';
v(8) = std(y(2:end)./y(1:end-1)-1)*sqrt(wid_year);%年化标准差
sta_val.std_return = v(8);
%9夏普比率：（策略日收益率的252日平均/策略日收益率的252日标准差）*[252^（1/2）]
%v(9) = v(1)/v(8)*sqrt(252);
v_str{9} = '夏普比率';
%v(9) = (((mean(y(2:end)./y(1:end-1)-1)-3/100/252))/std(y(2:end)./y(1:end-1)-1))*sqrt(252);
temp = y(2:end)./y(1:end-1)-1;
temp(isinf(temp)|isnan(temp)) = [];
%v(9) = ((mean(temp)-(exp(log(1.03)/252)-1)))/(std(temp))*sqrt(wid_year);
v(9) = (mean(temp))/(std(temp))*sqrt(wid_year); %20200504更新 无风险收益记为0
sta_val.sharp = v(9);
%10最长新高间隔
v_str{10} = '最长新高间隔';
temp = getdrawdown(y);
ind = find(eq(temp,0));
if length(ind)>1
    v(10) = max(diff(ind));
else
    v(10) = nan;
end
sta_val.max_inter = v(10);
%11 小于0的收益的波动率
v_str{11} = '下跌波动率';
v(11) = std(y(2:end)./y(1:end-1)-1);
sta_val.down_std = v(11);
v1 = y(wid_year:end)./y(1:end-wid_year+1)-1;
v_str{12} = '持有1年收益中位数';
v(12) = median(v1);
sta_val.median_return = v(12);
v_str{13} = '持有1年收益90%分位数';
v(13) = prctile(v1,90);
sta_val.ninety_quantile = v(13);

v_str{14} = '信息比率';
v(14) = v(1)/v(8);
sta_val.ninety_quantile = v(14);

temp = y(2:end)./y(1:end-1)-1;
temp(eq(temp,0)) = [];
v_str{15} = '胜率';
v(15) = sum(temp>0)/length(temp);
sta_val.ninety_quantile = v(14);


if disp_sel
str = [];
for i = 1:length(v)
    str = [str,sprintf('%s: %0.4f\n',v_str{i},v(i))];
end
sprintf('回测曲线参数：%s \n',str)
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
% %计算评价参数
% %年化收益-计算公式为：（账户最终价值/账户初始价值）^（250/回测期间总天数）-1
% value = [];
% value.nianhuashouyi = (MT(end)/MT(1))^(250/(N*Wnum*7))-1;
% %最大回撤-计算公式为：min(账户当日价值 / 当日之前账户最高价值-1)。
% tempv = zeros(size(MT));
% for i = 1:length(MT)
%     tempv(i) = MT(i)./max(MT(1:i))-1;
% end
% value.zuidahuiche = min(tempv);
% %平均涨幅 计算公式为：账户周收益的平均值。
% value.pingjunzhangfu = mean(MT(2:end)./MT(1:end-1)/Wnum-1)/Wnum;
% %上涨概率-计算公式为：上涨周数 / 回测交易周数量。
% value.shangzhanggailv = length(find(diff(MT)>0))/(length(MT)-1);
% %最大连续上涨周数 衡量了我们策略连续盈利的最大周数。
% testv1 = find(diff(MT)>0);
% %找出连续上涨的周序号
% testv2 = diff(testv1);
% %连续的序号差为1，否者大于1
% testv3 = [(1:size(testv2))',testv2];
% testv3(eq(testv2,1),:)=[];
% value.maxv1 = max(diff(testv3(:,1)));
% %最大连续下跌周数
% testv1 = find(diff(MT)<0);
% %找出下跌的周序号
% testv2 = diff(testv1);
% %连续的序号差为1，否者大于1
% testv3 = [(1:size(testv2))',testv2];
% testv3(eq(testv2,1),:)=[];
% value.maxv2 = max(diff(testv3(:,1)));
% %所有日期中最大的涨幅，衡量了资产一周内的最好表现。
% value.maxdanqizhangfu = max(MT(2:end)./MT(1:end-1)-1);
% %最大单周期跌幅
% value.maxdanqidiefu = min(MT(2:end)./MT(1:end-1)-1);
% %收益波动率
% value.shouyibodonglv = std(MT(2:end)./MT(1:end-1)-1)*sqrt(52/Wnum);%年化标准差
