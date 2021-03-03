function prt_sar(results,vnames,fid)
% PURPOSE: Prints output using SAR results structures
%---------------------------------------------------
% USAGE: prt_sar(results,vnames,fid)
% Where: results = a structure returned by a SAR model
%        vnames  = an optional vector of variable names
%        fid     = optional file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%--------------------------------------------------- 
%  NOTES: e.g. vnames = strvcat('y','const','x1','x2');
%         e.g. fid = fopen('ols.out','wr');
%  use prt_spat(results,[],fid) to print to a file with no vnames               
% --------------------------------------------------
%  RETURNS: nothing, just prints the SAR results
% --------------------------------------------------
% SEE ALSO: prt, plt
%---------------------------------------------------   

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

if ~isstruct(results)
 error('prt_sar requires structure argument');
elseif nargin == 1
 nflag = 0; fid = 1;
elseif nargin == 2
 fid = 1; nflag = 1;
elseif nargin == 3
 nflag = 0;
 [vsize junk] = size(vnames); % user may supply a blank argument
   if vsize > 0
   nflag = 1;          
   end;
else
 error('Wrong # of arguments to prt_sar');
end;


nvar = results.nvar;
nobs = results.nobs;

% handling of vnames
Vname = 'Variable';
 for i=1:nvar
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');

if (nflag == 1) % the user supplied variable names
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_sar -- check vnames argument \n');
 fprintf(fid,'will use generic variable names \n');
 nflag = 0;
 else,
Vname = 'Variable';
 for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
 end; % end of if-else
end; % end of nflag issue



switch results.meth

case {'sar'} % <=================== max lik spatial autoregressive model

nobs = results.nobs;
nvar = results.nvar;

fprintf(fid,'\n');
fprintf(fid,'Spatial autoregressive Model Estimates \n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared          = %9.4f \n',results.rsqr);
fprintf(fid,'Rbar-squared       = %9.4f \n',results.rbar);
fprintf(fid,'sigma^2            = %9.4f \n',results.sige);
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'log-likelihood     = %16.8g \n',results.lik);
fprintf(fid,'# of iterations    = %6d   \n',results.iter);
fprintf(fid,'min and max rho    = %9.4f,%9.4f \n',results.rmin,results.rmax);
% print timing information
fprintf(fid,'total time in secs = %9.4f \n',results.time);
if results.time1 ~= 0
fprintf(fid,'time for lndet     = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for eigs      = %9.4f \n',results.time2);
end;
if results.time3 ~= 0
fprintf(fid,'time for t-stats   = %9.4f \n',results.time3);
end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.miter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

fprintf(fid,'***************************************************************\n');

bout = [results.beta
        results.rho];
        
% <=================== end of sar case

case {'sar_g'} % <=================== MCMC spatial autoregressive model


nobs = results.nobs;
nvar = results.nvar;

% extract posterior means
bout = [results.beta
        results.rho];
sige = results.sige;
    tmp1 = std(results.bdraw);
    tmp2 = std(results.pdraw);
    bstd = [tmp1'
            tmp2];  

if strcmp(results.tflag,'tstat')
 tstat = bout./bstd;
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw];
 for i=1:results.nvar+1;
 if bout(i,1) > 0
 cnt = find(draws(:,i) > 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 else
 cnt = find(draws(:,i) < 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 end; % end of if - else
 end; % end of for loop
end; 

rsqr = results.rsqr;

fprintf(fid,'\n');
fprintf(fid,'Bayesian spatial autoregressive model \n');
if results.novi == 1
    fprintf(fid,'Homoscedastic version \n');
elseif results.novi == 0
    fprintf(fid,'Heteroscedastic model \n');
end;
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared          = %9.4f \n',rsqr);
fprintf(fid,'Rbar-squared       = %9.4f \n',results.rbar);
fprintf(fid,'mean of sige draws = %9.4f \n',results.sige);
fprintf(fid,'sige, epe/(n-k)    = %9.4f \n',results.sigma);
if (results.rdraw == 0 & results.novi == 0)
fprintf(fid,'r-value            = %6d   \n',results.r);
elseif (results.rdraw ~= 0  & results.novi == 0)
fprintf(fid,'mean of rdraws     = %9.4f \n',mean(results.rdraw));
fprintf(fid,'gam(m,k) prior     = %6d,%6d \n',results.m,results.k);
end;  
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'ndraws,nomit       = %6d,%6d \n',results.ndraw,results.nomit);
fprintf(fid,'total time in secs = %9.4f   \n',results.time);
if results.time1 ~= 0
fprintf(fid,'time for eigs      = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for lndet     = %9.4f \n',results.time2);
end;
if results.time3 ~= 0
fprintf(fid,'time for sampling  = %9.4f \n',results.time3);
end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.iter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

fprintf(fid,'min and max rho    = %9.4f,%9.4f \n',results.rmin,results.rmax);
fprintf(fid,'***************************************************************\n');


if (results.priorb == 1)
    % non-diffuse prior, so print it
vstring = 'Variable';
bstring = 'Prior Mean';
tstring = 'Std Deviation';

tmp = [results.bmean results.bstd];

cnames = strvcat(bstring,tstring);
rnames = vstring;
for i=1:nvar
rnames = strvcat(rnames,Vname(i+1,:));
end;

pin.fmt = '%16.6f';
pin.fid = fid;
pin.cnames = cnames;
pin.rnames = rnames;

mprint(tmp,pin);
fprintf(fid,'***************************************************************\n');
end;
fprintf(fid,'      Posterior Estimates \n');

 if strcmp(results.tflag,'tstat')
% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities
      
tmp = [bout results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
 else % use p-levels for Bayesian results
tmp = [bout bstd tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Std Deviation'; pstring = 'p-level';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
end;

return;

% <=================== end of sar_g case

case {'sar_gv'} % <=================== MCMC spatial autoregressive model
                % that produces r-value posterior based on draws


nobs = results.nobs;
nvar = results.nvar;

% extract posterior means
bout = [results.beta
        results.rho];
sige = results.sige;
    tmp1 = std(results.bdraw);
    tmp2 = std(results.pdraw);
    bstd = [tmp1'
            tmp2];  

if strcmp(results.tflag,'tstat')
 tstat = bout./bstd;
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw];
 for i=1:results.nvar+1;
 if bout(i,1) > 0
 cnt = find(draws(:,i) > 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 else
 cnt = find(draws(:,i) < 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 end; % end of if - else
 end; % end of for loop
end; 

rsqr = results.rsqr;

fprintf(fid,'\n');
fprintf(fid,'Bayesian spatial autoregressive model \n');
    fprintf(fid,'Heteroscedastic version with r-value estimate \n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared          = %9.4f \n',rsqr);
fprintf(fid,'Rbar-squared       = %9.4f \n',results.rbar);
fprintf(fid,'mean of sige draws = %9.4f \n',results.sige);
fprintf(fid,'sige, epe/(n-k)    = %9.4f \n',results.sigma);
fprintf(fid,'mean r-value       = %9.4f \n',mean(results.rdraw));
fprintf(fid,'std  r-value       = %9.4f \n',std(results.rdraw));
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'ndraws,nomit       = %6d,%6d \n',results.ndraw,results.nomit);
fprintf(fid,'total time in secs = %9.4f   \n',results.time);
if results.time1 ~= 0
fprintf(fid,'time for eigs      = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for lndet     = %9.4f \n',results.time2);
end;
if results.time3 ~= 0
fprintf(fid,'time for sampling  = %9.4f \n',results.time3);
end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.iter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

fprintf(fid,'min and max rho    = %9.4f,%9.4f \n',results.rmin,results.rmax);
fprintf(fid,'***************************************************************\n');


if (results.priorb == 1)
    % non-diffuse prior, so print it
vstring = 'Variable';
bstring = 'Prior Mean';
tstring = 'Std Deviation';

tmp = [results.bmean results.bstd];

cnames = strvcat(bstring,tstring);
rnames = vstring;
for i=1:nvar
rnames = strvcat(rnames,Vname(i+1,:));
end;

pin.fmt = '%16.6f';
pin.fid = fid;
pin.cnames = cnames;
pin.rnames = rnames;

mprint(tmp,pin);
fprintf(fid,'***************************************************************\n');
end;
fprintf(fid,'      Posterior Estimates \n');

 if strcmp(results.tflag,'tstat')
% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities
      
tmp = [bout results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
 else % use p-levels for Bayesian results
tmp = [bout bstd tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Std Deviation'; pstring = 'p-level';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
end;

return;

% <=================== end of sar_gv case


case {'sar_c'} % <=================== log-marginal for spatial autoregressive model

fprintf(fid,'sar_c: no printed output available, this function just produces a log-marginal estimates \n');

% ,============ end of sar_c case

case {'sart_g','sart_gc'} % <=================== Gibbs spatial autoregressive Tobit model

nobs = results.nobs;
nvar = results.nvar;


% find posterior means
tmp1 = mean(results.bdraw);
pout = mean(results.pdraw);
bout = [tmp1'
        pout];

y = results.y;
yhat = results.yhat;
sige = mean(results.sdraw);
tmp1 = std(results.bdraw);
tmp2 = std(results.pdraw);
bstd = [tmp1'
        tmp2];

if strcmp(results.tflag,'tstat')
 tstat = bout./bstd;
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw];
 for i=1:results.nvar+1;
 if bout(i,1) > 0
 cnt = find(draws(:,i) > 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 else
 cnt = find(draws(:,i) < 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 end; % end of if - else
 end; % end of for loop
end; 



fprintf(fid,'\n');
fprintf(fid,'Bayesian spatial autoregressive Tobit model \n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'mean of sige draws = %9.4f \n',sige);
if results.rdraw == 0
fprintf(fid,'r-value            = %6d   \n',results.r);
else
fprintf(fid,'mean of rdraws     = %9.4f \n',mean(results.rdraw));
fprintf(fid,'gam(m,k) prior     = %6d,%6d \n',results.m,results.k);
end; 
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'# censored values  = %6d \n',results.nobsc);
fprintf(fid,'ndraws,nomit       = %6d,%6d \n',results.ndraw,results.nomit);
fprintf(fid,'time in secs       = %9.4f   \n',results.time);
fprintf(fid,'min and max rho    = %9.4f,%9.4f \n',results.rmin,results.rmax);
fprintf(fid,'***************************************************************\n');

if (results.priorb == 1)
    % non-diffuse prior, so print it
vstring = 'Variable';
bstring = 'Prior Mean';
tstring = 'Std Deviation';

tmp = [results.bmean results.bstd];

cnames = strvcat(bstring,tstring);
rnames = vstring;
for i=1:nvar
rnames = strvcat(rnames,Vname(i+1,:));
end;

pin.fmt = '%16.6f';
pin.fid = fid;
pin.cnames = cnames;
pin.rnames = rnames;

mprint(tmp,pin);
fprintf(fid,'***************************************************************\n');
end;

fprintf(fid,'      Posterior Estimates \n');

 if strcmp(results.tflag,'tstat')
% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities
      
tmp = [bout results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
 else % use p-levels for Bayesian results
tmp = [bout bstd tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Std Deviation'; pstring = 'p-level';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
end;

return;

% <=================== end of sart_g case


case {'sarp_g','sarp_gc'} % <=================== Gibbs spatial autoregressive Probit model

nobs = results.nobs;
nvar = results.nvar;
  
% find posterior means
tmp1 = mean(results.bdraw);
pout = mean(results.pdraw);
bout = [tmp1'
        pout];

tmp1 = std(results.bdraw);
tmp2 = std(results.pdraw);
bstd = [tmp1'
        tmp2];

if strcmp(results.tflag,'tstat')
 tstat = bout./bstd;
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw];
 for i=1:results.nvar+1;
 if bout(i,1) > 0
 cnt = find(draws(:,i) > 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 else
 cnt = find(draws(:,i) < 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 end; % end of if - else
 end; % end of for loop
end; 


fprintf(fid,'\n');
fprintf(fid,'Bayesian spatial autoregressive Probit model \n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'psuedo R-sqr    = %9.4f \n',results.rsqr);
fprintf(fid,'sige            = %9.4f \n',mean(results.sdraw));
fprintf(fid,'Nobs, Nvars     = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'# 0, 1 y-values = %6d,%6d \n',results.zip,nobs-results.zip);
fprintf(fid,'ndraws,nomit    = %6d,%6d \n',results.ndraw,results.nomit);
fprintf(fid,'time in secs    = %9.4f   \n',results.time);
fprintf(fid,'min and max rho = %9.4f,%9.4f \n',results.rmin,results.rmax);
fprintf(fid,'***************************************************************\n');


 if strcmp(results.tflag,'tstat')
% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities
      
tmp = [bout results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
 else % use p-levels for Bayesian results
tmp = [bout bstd tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Std Deviation'; pstring = 'p-level';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
end;

return;

% <=================== end of sarp_g case

case {'sar_gbma'} % <=================== bma for sar models
nmodels = results.nmodels;
nvar = results.nvar;
if nargin < 3
fid = 1;
end;

mout = results.modelsa;
occ = sum(mout);

fmt = [];
for i=1:nvar;
fmt = strvcat(fmt,'%5d');
end;
fmt = strvcat(fmt,'%8.4f');
in.fmt = fmt;
% 
rnames = 'Model';
for i=1:nmodels;
rnames = strvcat(rnames,['model ' num2str(i)]); 
end;
rnames = strvcat(rnames,'#Occurences');
in.rnames = rnames;
cnames = vnames(2:end,:);
cnames = strvcat(cnames,'probs');
in.cnames = cnames;
% 
fprintf(fid,'Model averaging information \n');
in.width = 3000;
in.fid = fid;
[tst1,tst2] = size(results.models);
if tst1 > 1 
out = [results.models
       occ];
mprint(out,in);
else
out = [results.models];
in.rnames = rnames(end-1:end,:);
end;
mprint(out,in);
fprintf(fid,'***************************************************************\n');

% only do this is avg_flag == 1
if results.avg_flag == 1
sige = mean(results.sdraw);
bhat = mean(results.bdraw);
bstd = std(results.bdraw);
bhatp = bhat';
bstdp = bstd';
rho = mean(results.pdraw);
rstd = std(results.pdraw);

fprintf(fid,'\n');
fprintf(fid,'SAR Bayesian Model Averaging Estimates \n');
fprintf(fid,'Dependent Variable   = %16s \n',vnames(1,:));
fprintf(fid,'R-squared            = %9.4f \n',results.rsqr);
fprintf(fid,'sigma^2              = %9.4f \n',sige);
fprintf(fid,'# unique models      = %10d \n',results.munique);
fprintf(fid,'Nobs, Nvars          = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'ndraws for BMA       = %6d \n',results.ndraw);
fprintf(fid,'ndraws for estimates = %6d \n',results.ndraw2);
fprintf(fid,'nomit for estimates  = %6d \n',results.nomit2);

if results.time1 ~= 0
fprintf(fid,'time for eigs        = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for lndet       = %9.4f \n',results.time2);
end;
if results.time3 ~= 0
fprintf(fid,'time for BMA sampling= %9.4f \n',results.time3);
end;
if results.time4 ~= 0
fprintf(fid,'time for estimates   = %9.4f \n',results.time4);
end;


if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.iter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

fprintf(fid,'min and max rho    = %9.4f,%9.4f \n',results.rmin,results.rmax);

vstring = 'Variable';
bstring = 'Prior Mean';
tstring = 'Std Deviation';

tmp = [results.bmean results.bstd];

cnames = strvcat(bstring,tstring);
rnames = vstring;
for i=1:nvar
rnames = strvcat(rnames,Vname(i+1,:));
end;

pin.fmt = '%16.6f';
pin.fid = fid;
pin.cnames = cnames;
pin.rnames = rnames;
fprintf(fid,'***************************************************************\n');

mprint(tmp,pin);

fprintf(fid,'***************************************************************\n');
fprintf(fid,'      Posterior Estimates \n');


% column labels for printing results
vstring = 'Variable';
bstring = strvcat('Coefficient','std dev');

ball = [bhatp bstdp
        rho   rstd];


if strcmp(results.tflag,'tstat')
% now print coefficient estimates, t-statistics and probabilities
tstat = [bhatp./bstdp
         rho/rstd];
tout = norm_prb(tstat); % find asymptotic z (normal) probabilities    
tmp = [ball(:,1) tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
rnames = 'Variable';
 for i=1:nvar
  rnames = strvcat(rnames,vnames(i+1,:));
 end;
rnames = strvcat(rnames,'rho');

in.cnames = cnames;
in.rnames = rnames;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
 else % use p-levels for Bayesian results
 draws = [results.bdraw results.pdraw];

 for i=1:results.nvar+1;
 if ball(i,1) > 0
 cnt = find(draws(:,i) > 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw2 - results.nomit2));
 else
 cnt = find(draws(:,i) < 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw2 - results.nomit2));
 end; % end of if - else
 end; % end of for loop

tmp = [ball tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Std Deviation'; pstring = 'p-level';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
rnames = 'Variable';
 for i=1:nvar
  rnames = strvcat(rnames,vnames(i+1,:));
 end;
rnames = strvcat(rnames,'rho');
in.rnames = rnames;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
end;

end; % end of avg_flag == 1

return;


otherwise
error('results structure not known by prt_sar function');
end;


% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities
tmp = [bout results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
