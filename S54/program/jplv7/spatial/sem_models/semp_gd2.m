% PURPOSE: An example of using semp_g()
%          Gibbs sampling spatial autoregressive model
%          on a large data set                    
%---------------------------------------------------
% USAGE: semp_gd2 (see semp_gd for a small data set)
%---------------------------------------------------

clear all;

load election.data;
%   documentation for election.data
%   3,110 Counties
%   
%   education variables are expressed as a proportion of all those with degrees,
%   that is as a proportion of: educated = (high+smcollege+associate+college+gradprof)
%   NOTE:because of this, the education proportion variables sum to unity
%   
%   col 1= binary y with 0=Dole, 1=Clinton (1996 Presidential Election)
%   col2 = latt  coordinate
%   col3 = long  coordinate
%   col4 = log-urban  (log of urban population)      
%   col5 = log-rural  (log of rural population)    
%   col6 = prop-highs (pop with high school or GED graduates as aproportion of educated -- see above )  
%   col7 = prop-smcollege  (pop with some college as aproportion of educated -- see above)  
%   col8 = prop-associate (pop with associate degrees as aproportion of educated -- see above)      
%   col9 = prop-college   (pop with college degrees as aproportion of educated -- see above)  
%   col10 = prop-gradprof (pop with grad/professional degrees as aproportion of educated -- see above) 
%   col11 = statecode (a number 1 to 48 with the state in which the county is located)

y = election(:,1);
xc = election(:,2);
yc = election(:,3);
n = length(y);

xmat = [ones(n,1) election(:,4) election(:,7:10)];

[j1 W j2] = xy2cont(xc,yc);
clear election;                % conserve on RAM memory
n = 3110;
vnames = strvcat('clinton/dole,s','constant','urban','some college', ...
'associate degree','college','graduate/professional');

ndraw = 2500; 
nomit = 500;

result0 = probit(y,xmat);
prt(result0,vnames);

prior.rval = 4;
result = semp_g(y,xmat,W,ndraw,nomit,prior);
prt(result,vnames);

tt=1:n;
plot(tt,y,'og',tt,result.yprob,'.r',tt,result0.yhat,'.b');
legend('0,1 y-values','sem probs','probit probs');
pause;

prior2.dflag = 1;
prior2.rval = 4;
result2 = sarp_g(y,xmat,W,ndraw,nomit,prior2);
prt(result2,vnames);

tt=1:n;
plot(tt,y,'og',tt,result2.yprob,'.r');
legend('0,1 y-values','sar probs');

% compute posterior model probabilities
probs = model_probs(result,result2);

fprintf(1,'posterior probs for sarp versus semp model \n');
in.rnames = strvcat('Models','sem_probit','sar_probit');
mprint(probs,in);

