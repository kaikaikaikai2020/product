%将xlsx转换为csv格式
function xlstocsv_adair(fn,x)

[pn,sub_fn] = fileparts(fn);
app_r = '.csv';
fn = fullfile(pn,[sub_fn,app_r]);

if iscell(x)
    x = cell2table(x);
elseif isnumeric(x)
    x = array2table(x);
end

writetable(x,fn)