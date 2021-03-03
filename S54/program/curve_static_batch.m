 function sta_re = curve_static_batch(yc,title_str)
    if ~iscell(yc)
        temp = cell(size(yc(1,:)));
        for i = 1:size(yc,2)
            temp{i} = yc(:,i);
        end
        yc = temp;
    end
    sta_re = cell(size(yc));
    for i = 1:length(sta_re)
        [v0,v_str0] = curve_static(yc{i},[]);
        [v,v_str] = ad_trans_sta_info(v0,v_str0);
        if eq(i,1)
            sub_re = [[{''};v_str'],[title_str(i);v']];
        else
            sub_re = [title_str(i);v'];
        end
        sta_re{i} = sub_re;
    end
    sta_re = [sta_re{:}]';
end