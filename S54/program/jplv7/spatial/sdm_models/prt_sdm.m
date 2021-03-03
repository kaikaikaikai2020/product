function prt_sdm(results,vnames,fid)
% PURPOSE: Prints output using sdm results structures
%---------------------------------------------------
% USAGE: prt_sdm(results,vnames,fid)
% Where: results = a structure returned by sdm, sdm_g, sdm_gc, etc.
%        vnames  = an optional vector of variable names
%        fid     = optional file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%--------------------------------------------------- 
%  NOTES: e.g. vnames = strvcat('y','const','x1','x2');
%         e.g. fid = fopen('ols.out','wr');
%  use prt_sdm(results,[],fid) to print to a file with no vnames               
% --------------------------------------------------
%  RETURNS: nothing, just prints the spatial regression results
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
 error('prt_sdm requires structure argument');
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
 error('Wrong # of arguments to prt_sdm');
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
 fprintf(fid,'Wrong # of variable names in prt_sdm -- check vnames argument \n');
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

case {'sdm'} % <=================== spatial durbin model

nobs = results.nobs;
nvar = results.nvar;

% special handling of vnames
Vname = 'Variable';
if results.cflag == 1 % we have an intercept
 for i=1:2*(nvar-1)+1;
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
elseif results.cflag == 0 % no intercept term
 for i=1:2*(nvar);
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
end;

if (nflag == 1) % the user supplied variable names
Vname = 'Variable';
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_sdm -- check vnames argument \n');
 nflag = 0;
 fprintf(fid,'will use generic variable names \n');
 else
  if results.cflag == 1 % we have an intercept
   for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
   end;
   for i=2:nvar
    Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
   end;
   % add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
   elseif results.cflag == 0
   for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
   end;
   for i=1:nvar
    Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
   end;
   % add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
   end; % end of if/elseif
  end; % end of if-else
end; % end of nflag issue


fprintf(fid,'\n');
fprintf(fid,'Spatial Durbin model\n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared          = %9.4f   \n',results.rsqr);
fprintf(fid,'Rbar-squared       = %9.4f   \n',results.rbar);
fprintf(fid,'sigma^2            = %9.4f   \n',results.sige);
fprintf(fid,'log-likelihood     = %16.8g  \n',results.lik);
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'# iterations       = %6d     \n',results.iter);
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
        
% <=================== end of sdm case

case {'sdm_g'} % <=================== spatial durbin model MCMC


nobs = results.nobs;
nvar = results.nvar;

Vname = 'Variable';
if results.cflag == 1 % we have an intercept
 for i=1:2*(nvar-1)+1;
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
elseif results.cflag == 0 % no intercept term
 for i=1:2*(nvar);
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
end;

if (nflag == 1) % the user supplied variable names
Vname = 'Variable';
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_sdm -- check vnames argument \n');
 nflag = 0;
 fprintf(fid,'will use generic variable names \n');
 else
  if results.cflag == 1 % we have an intercept
   for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
   end;
   for i=2:nvar
    Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
   end;
   % add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
   elseif results.cflag == 0
   for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
   end;
   for i=1:nvar
    Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
   end;
   % add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
   end; % end of if/elseif
  end; % end of if-else
end; % end of nflag issue


% find posterior means
    bout = [results.beta
            results.rho];
    sige = results.sige;
    tmp1 = std(results.bdraw);
    tmp2 = std(results.pdraw);
    bstd = [tmp1'
        tmp2];  


if strcmp(results.tflag,'tstat')
 tstat = bout./bstd;
 [junk nk] = size(results.bdraw);
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw];
 [junk nk] = size(draws);
 for i=1:nk;
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
fprintf(fid,'Bayesian Spatial Durbin model\n');
if results.novi == 1
    fprintf(fid,'Homoscedastic version \n');
elseif results.novi == 0
    fprintf(fid,'Heteroscedastic model \n');
end;    
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared          = %9.4f   \n',results.rsqr);
fprintf(fid,'mean of sige draws = %9.4f   \n',sige);
fprintf(fid,'sige, epe/(n-k)    = %9.4f   \n',results.sigma);
if (results.rdraw == 0 & results.novi == 0)
fprintf(fid,'r-value            = %6d   \n',results.r);
elseif (results.rdraw ~= 0  & results.novi == 0)
fprintf(fid,'mean of rdraws     = %9.4f \n',mean(results.rdraw));
fprintf(fid,'gam(m,k) prior     = %6d,%6d \n',results.m,results.k);
end;  
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,nk);
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
fprintf(fid,'min and max rho= %9.4f,%9.4f \n',results.rmin,results.rmax);
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
   
% <=================== end of sdm_g case
     
case {'sdm_gc'} % <=================== spatial durbin model MCMC

nobs = results.nobs;
nvar = results.nvar;

Vname = 'Variable';
if results.cflag == 1 % we have an intercept
 for i=1:2*(nvar-1)+1;
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
elseif results.cflag == 0 % no intercept term
 for i=1:2*(nvar);
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
end;

if (nflag == 1) % the user supplied variable names
Vname = 'Variable';
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_sdm -- check vnames argument \n');
 nflag = 0;
 fprintf(fid,'will use generic variable names \n');
 else
  for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
  end;
  for i=2:nvar
    Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
  end;
  % add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
 end; % end of if-else
end; % end of nflag issue


% find posterior means
    tmp1 = mean(results.bdraw);
    pout = mean(results.pdraw);
    bout = [tmp1'
        pout];
    sige = mean(results.sdraw);
    tmp1 = std(results.bdraw);
    tmp2 = std(results.pdraw);
    bstd = [tmp1'
        tmp2];  


if strcmp(results.tflag,'tstat')
 tstat = bout./bstd;
 [junk nk] = size(results.bdraw);
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw];
 [junk nk] = size(draws);
 for i=1:nk;
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
fprintf(fid,'Bayesian Spatial Durbin model\n');
if results.novi == 1
    fprintf(fid,'Homoscedastic version \n');
elseif results.novi == 0
    fprintf(fid,'Heteroscedastic model \n');
end;
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared          = %9.4f   \n',results.rsqr);
fprintf(fid,'sigma^2            = %9.4f   \n',results.sige);
if (results.rdraw == 0 & results.novi == 0)
fprintf(fid,'r-value            = %6d   \n',results.r);
elseif (results.rdraw ~= 0  & results.novi == 0)
fprintf(fid,'mean of rdraws     = %9.4f \n',mean(results.rdraw));
fprintf(fid,'gam(m,k) prior     = %6d,%6d \n',results.m,results.k);
end;  
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,nk);
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

fprintf(fid,'min and max rho= %9.4f,%9.4f \n',results.rmin,results.rmax);
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
   
% <=================== end of sdm_gc case

case {'sdmp_g','sdmp_gc'} % <=================== spatial durbin probit model MCMC

nobs = results.nobs;
nvar = results.nvar;

Vname = 'Variable';
if results.cflag == 1 % we have an intercept
 for i=1:2*(nvar-1)+1;
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
elseif results.cflag == 0 % no intercept term
 for i=1:2*(nvar);
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
end;

if (nflag == 1) % the user supplied variable names
Vname = 'Variable';
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_sdm -- check vnames argument \n');
 nflag = 0;
 fprintf(fid,'will use generic variable names \n');
 else
  if results.cflag == 1 % we have an intercept
   for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
   end;
   for i=2:nvar
    Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
   end;
   % add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
   elseif results.cflag == 0
   for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
   end;
   for i=1:nvar
    Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
   end;
   % add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
   end; % end of if/elseif
  end; % end of if-else
end; % end of nflag issue


% find posterior means
    tmp1 = mean(results.bdraw);
    pout = mean(results.pdraw);
    bout = [tmp1'
        pout];
    sige = mean(results.sdraw);
    tmp1 = std(results.bdraw);
    tmp2 = std(results.pdraw);
    bstd = [tmp1'
        tmp2];  


if strcmp(results.tflag,'tstat')
 tstat = bout./bstd;
 [junk nk] = size(results.bdraw);
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw];
 [junk nk] = size(draws);
 for i=1:nk;
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
fprintf(fid,'Bayesian Spatial Durbin Probit model\n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'sigma^2            = %9.4f   \n',sige);
fprintf(fid,'# 0, 1 y-values    = %6d,%6d \n',results.zip,nobs-results.zip);

if results.rdraw == 0
fprintf(fid,'r-value            = %6d   \n',results.r);
else
fprintf(fid,'mean of rdraws     = %9.4f \n',mean(results.rdraw));
fprintf(fid,'gam(m,k) prior     = %6d,%6d \n',results.m,results.k);
end;    
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,nk);
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

fprintf(fid,'min and max rho= %9.4f,%9.4f \n',results.rmin,results.rmax);
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
   
% <=================== end of sdmp_g case
     
case {'sdmt_g','sdmt_gc'} % <=================== spatial durbin tobit model MCMC

nobs = results.nobs;
nvar = results.nvar;

Vname = 'Variable';
if results.cflag == 1 % we have an intercept
 for i=1:2*(nvar-1)+1;
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
elseif results.cflag == 0 % no intercept term
 for i=1:2*(nvar);
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
end;

if (nflag == 1) % the user supplied variable names
Vname = 'Variable';
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_sdm -- check vnames argument \n');
 nflag = 0;
 fprintf(fid,'will use generic variable names \n');
 else
  if results.cflag == 1 % we have an intercept
   for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
   end;
   for i=2:nvar
    Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
   end;
   % add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
   elseif results.cflag == 0
   for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
   end;
   for i=1:nvar
    Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
   end;
   % add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
   end; % end of if/elseif
  end; % end of if-else
end; % end of nflag issue


% find posterior means
    tmp1 = mean(results.bdraw);
    pout = mean(results.pdraw);
    bout = [tmp1'
        pout];
    sige = mean(results.sdraw);
    tmp1 = std(results.bdraw);
    tmp2 = std(results.pdraw);
    bstd = [tmp1'
        tmp2];  


if strcmp(results.tflag,'tstat')
    [junk nk] = size(results.bdraw);
 tstat = bout./bstd;
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw];
 [junk nk] = size(draws);
 for i=1:nk;
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
fprintf(fid,'Bayesian Spatial Durbin tobit model\n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'sige               = %9.4f \n',sige);
if results.rdraw == 0
fprintf(fid,'r-value            = %6d   \n',results.r);
else
fprintf(fid,'mean of rdraws     = %9.4f \n',mean(results.rdraw));
fprintf(fid,'gam(m,k) prior     = %6d,%6d \n',results.m,results.k);
end;    
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,nk);
fprintf(fid,'# censored values  = %6d \n',results.nobsc);
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

fprintf(fid,'min and max rho= %9.4f,%9.4f \n',results.rmin,results.rmax);
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
   
% <=================== end of sdmt_g case

otherwise
error('results structure not known by prt_sdm function');
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
