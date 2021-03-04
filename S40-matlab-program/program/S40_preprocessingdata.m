function  [t,testStart,openprice,closeprice] = S40_preprocessingdata(x,N)
if nargin < 2
    N = 2;
end
temp_num = x(:,4)*100+x(:,5);
id = temp_num>=930&temp_num<=1500;
x = x(id,:);
t = datenum([x(:,1:5),zeros(size(x(:,1)))]);
%统计中间有停牌的情况，并剔除
day_tick = x(:,1)*10000+x(:,2)*100+x(:,3);
day_tick_u = unique(day_tick);
ind_miss = false(size(day_tick));
ind_miss_u = false(size(day_tick_u));
T = length(day_tick_u);
for i = 1:T
    sub_ind = eq(day_tick,day_tick_u(i));
    if sum(sub_ind)<240
        ind_miss(sub_ind) = true;
        ind_miss_u(i) = true;
    end
end

day_tick_u(ind_miss_u) = [];
x(ind_miss,:) = [];
t(ind_miss,:) = [];
%初始时间设定
min_day_num = 210*N;
temp = num2str(day_tick_u(min_day_num/2+1));
testStart=sprintf('%s-%s-%s',temp(1:4),temp(5:6),temp(7:8));

openprice = x(:,end-1);
closeprice = x(:,end);
end