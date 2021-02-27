function r = get_ret(sub_pos1,sub_pos2,sub_fee_com)
    if nargin < 3
        sub_fee_com = zeros(length(sub_pos1),2);
    end
    r = zeros(1,2);
    for i = 1:2
        sub_p1 = sub_pos1(i);
        sub_p2 = sub_pos2(i);
        sub_fee = sub_fee_com(:,i);
        r(i) = get_ret0(sub_p1,sub_p2,sub_fee);
    end
end

