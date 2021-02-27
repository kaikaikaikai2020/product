function r = get_ret0(sub_p1,sub_p2,fees)
    if nargin < 3
        fees = [0,0];
    end
    
    fee_buy = fees(1);
    fee_sail = fees(2);
    
    if sub_p1>0
%         sub_p1 = sub_p1*(1+fee_buy);
%         sub_p2 = sub_p2*(1-fee_sail);
        r = (sub_p2-sub_p1)/sub_p1;
        r = r-fee_buy -fee_sail;
    elseif sub_p1<0
        %sub_p1 = sub_p1*(1-fee_sail);
        r = (abs(sub_p1)-sub_p2)/abs(sub_p1);
        r = r-fee_buy-fee_sail;
    else
        r = 0;
    end
end