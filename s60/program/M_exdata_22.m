clear

addpath(genpath(fullfile(pwd,'jplv7')))

load tempdata20200920
%usdtwd usdcnh
symbol_pool = {'NTN_1M','KWN_1M','NTN_6M','NTN_12M'};
T = length(symbol_pool);
for i = 1:T
    fn = symbol_pool{i};
    x = data.(fn);
    x = cell2table(x,'VariableNames',{'date','close'});
    writetable(x,sprintf('%s.csv',fn));
end