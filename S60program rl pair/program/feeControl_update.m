function [feeH,feeA,stampH,stampA,sRateH,sRateA]=feeControl_update(fee)
%% 函数说明
% 控制费率函数,ee=1考虑费用,fee=0不考虑费用
%%
if fee==1
    feeH=1/10000/2;
    feeA=1/10000/2;
    stampH=1/10000/2;
    stampA=1/10000/2;
    %sRateH=3/100;
    %sRateA=8.6/100;
    sRateH = 0;
    sRateA = 0;
else
    feeH=0/1000;
    feeA=0/1000;
    stampH=0/1000;
    stampA=0/1000;
    sRateH=0/100;
    sRateA=0/100;
end
end