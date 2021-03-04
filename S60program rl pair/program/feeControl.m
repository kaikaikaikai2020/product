function [feeH,feeA,stampH,stampA,sRateH,sRateA]=feeControl(fee,cal_mod)
if nargin < 2
    cal_mod=1;
end
%% 函数说明
% 控制费率函数,ee=1考虑费用,fee=0不考虑费用
%%
if fee==1
    feeH=1/1000;
    feeA=1/1000;
    stampH=1/1000;
    stampA=1/1000;
    if eq(cal_mod,1)
        sRateH=3/100;
        sRateA=8.6/100;
    else
        sRateH = 0;
        sRateA = 0;
    end
else
    feeH=0/1000;
    feeA=0/1000;
    stampH=0/1000;
    stampA=0/1000;
    sRateH=0/100;
    sRateA=0/100;
end
end