function [v,v_str0] = ad_trans_sta_info(v0,v_str0)
    v=v0;
    id1 = [1:5,8,11,12,13,15];
    v(id1) = v0(id1) * 100;
    v = num2cell(v);
    v_str0(id1) = cellfun(@(x) [x,'(%)'],v_str0(id1),'UniformOutput',false);
    id2 = [6,7,10];
    v(id2) = cellfun(@(x) sprintf('%d',x),v(id2),'UniformOutput',false);
    id3 = 1:length(v);
    id3(id2) = [];
    v(id3) = cellfun(@(x) sprintf('%0.2f',x),v(id3),'UniformOutput',false);
end