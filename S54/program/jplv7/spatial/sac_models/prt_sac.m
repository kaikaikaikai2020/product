function prt_sac(results,vnames,fid)
% PURPOSE: Prints output using sac model results structures
%---------------------------------------------------
% USAGE: prt_sac(results,vnames,fid)
% Where: results = a structure returned by sac, sac_g, etc.
%        vnames  = an optional vector of variable names
%        fid     = optional file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%--------------------------------------------------- 
%  NOTES: e.g. vnames = strvcat('y','const','x1','x2');
%         e.g. fid = fopen('ols.out','wr');
%  use prt_spat(results,[],fid) to print to a file with no vnames               
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
% jpl@jpl.econ.utoledo.edu

if ~isstruct(results)
 error('prt_sac requires structure argument');
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
 error('Wrong # of arguments to prt_sac');
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
 fprintf(fid,'Wrong # of variable names in prt_sac -- check vnames argument \n');
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


case {'sac'} % <=================== general spatial model

nobs = results.nobs;
nvar = results.nvar;

% special handling of vnames
if ( nflag == 0) %  no variable names supplied, make some up
Vname = 'Variable';
 for i=1:nvar
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho, lambda parameter name
    Vname = strvcat(Vname,'rho');
    Vname = strvcat(Vname,'lambda');

elseif (nflag == 1) % the user supplied variable names
Vname = 'Variable';
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 error('Wrong # of variable names in prt_spat -- check vnames argument');
 end;
 for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
 end;
% add spatial rho, lambda parameter name
    Vname = strvcat(Vname,'rho');
    Vname = strvcat(Vname,'lambda');

end; % end of nflag issue


fprintf(fid,'\n');
fprintf(fid,'General Spatial Model Estimates \n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared          = %9.4f \n',results.rsqr);
fprintf(fid,'Rbar-squared       = %9.4f \n',results.rbar);
fprintf(fid,'sigma^2            = %9.4f \n',results.sige);
fprintf(fid,'log-likelihood     = %16.8g \n',results.lik);
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'# iterations       = %6d \n',results.iter);
% print timing information
fprintf(fid,'total time in secs = %9.4f \n',results.time);
fprintf(fid,'time for optimiz   = %9.4f \n',results.time4);
if results.time1 ~= 0
fprintf(fid,'time for lndet     = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for eigs      = %9.4f \n',results.time2);
end;
if results.time3 ~= 0
fprintf(fid,'time for t-stat    = %9.4f \n',results.time3);
end;
fprintf(fid,'***************************************************************\n');

bout = [results.beta
        results.rho
        results.lam];
        

% <=================== end of sac case

case {'sac_g'} % <=================== Gibbs general spatial model

nobs = results.nobs;
nvar = results.nvar;

% handling of vnames
Vname = 'Variable';
 for i=1:nvar
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
% add spatial rho parameter name
    Vname = strvcat(Vname,'lambda');

if (nflag == 1) % the user supplied variable names
Vname = 'Variable';
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_spat -- check vnames argument \n');
 fprintf(fid,'will use generic variable names \n');
 nflag = 0;
 else,
 for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    Vname = strvcat(Vname,'lambda');
 end; % end of if-else
end; % end of nflag issue

  
% pull out posterior means
bmean = results.beta;
rho = results.rho;
lam = results.lam;
bout = [bmean
        rho
        lam];

y = results.y;
sige = results.sige;
tmp1 = std(results.bdraw);
tmp2 = std(results.pdraw);
tmp3 = std(results.ldraw);
bstd = [tmp1'
        tmp2
        tmp3];

if strcmp(results.tflag,'tstat')
 tstat = bout./bstd;
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw results.ldraw];
 for i=1:results.nvar+2;
 if bout(i,1) > 0
 cnt = find(draws(:,i) > 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 else
 cnt = find(draws(:,i) < 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 end; % end of if - else
 end; % end of for loop
end; 

e = y - results.yhat;
sigu = e'*e;
ym = y - mean(y);
rsqr1 = sigu;
rsqr2 = ym'*ym;
rsqr = 1.0 - rsqr1/rsqr2; % conventional r-squared

fprintf(fid,'\n');

fprintf(fid,'Bayesian general spatial model \n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared           = %9.4f \n',rsqr);
fprintf(fid,'sigma^2             = %9.4f \n',sige);
fprintf(fid,'Nobs, Nvars         = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'ndraws,nomit        = %6d,%6d \n',results.ndraw,results.nomit);
fprintf(fid,'time in secs        = %9.4f   \n',results.time);
fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
fprintf(fid,'min and max lambda  = %9.4f,%9.4f \n',results.lmin,results.lmax);
if (results.rdraw == 0)
fprintf(fid,'r-value             = %6d   \n',results.r);
elseif (results.rdraw ~= 0)
fprintf(fid,'mean of rdraws      = %9.4f \n',mean(results.rdraw));
fprintf(fid,'gam(m,k) prior      = %6d,%6d \n',results.m,results.k);
end;  

if (results.bflag == 1)
fprintf(fid,'***************************************************************\n');
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
end;

fprintf(fid,'***************************************************************\n');

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

% end of sac_g case

case {'sact_g'} % <=================== Gibbs general spatial Tobit model

nobs = results.nobs;
nvar = results.nvar;

% handling of vnames
Vname = 'Variable';
 for i=1:nvar
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
% add spatial rho parameter name
    Vname = strvcat(Vname,'lambda');

if (nflag == 1) % the user supplied variable names
Vname = 'Variable';
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_spat -- check vnames argument \n');
 fprintf(fid,'will use generic variable names \n');
 nflag = 0;
 else,
 for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    Vname = strvcat(Vname,'lambda');
 end; % end of if-else
end; % end of nflag issue

  
% find posterior means
tmp1 = mean(results.bdraw);
pout = mean(results.pdraw);
lout = mean(results.ldraw);
bout = [tmp1'
        pout
        lout];

y = results.y;
sige = mean(results.sdraw);
tmp1 = std(results.bdraw);
tmp2 = std(results.pdraw);
tmp3 = std(results.ldraw);
bstd = [tmp1'
        tmp2
        tmp3];

if strcmp(results.pflag,'tstat')
 tstat = bout./bstd;
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw results.ldraw];
 for i=1:results.nvar+2;
 if bout(i,1) > 0
 cnt = find(draws(:,i) > 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 else
 cnt = find(draws(:,i) < 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 end; % end of if - else
 end; % end of for loop
end; 

yhat = results.yhat;
e = y - yhat;
sigu = e'*e;
ym = y - mean(y);
rsqr1 = sigu;
rsqr2 = ym'*ym;
rsqr = 1.0 - rsqr1/rsqr2; % conventional r-squared
e = results.ymean - yhat;
sigu = e'*e;
ym = results.ymean - mean(results.ymean);
rsqr1 = sigu;
rsqr2 = ym'*ym;
rsqri = 1.0 - rsqr1/rsqr2; % non-conventional r-squared

fprintf(fid,'\n');
fprintf(fid,'Bayesian general spatial Tobit model \n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared          = %9.4f \n',rsqr);
fprintf(fid,'Imputed R-squared  = %9.4f \n',rsqri);
fprintf(fid,'sigma^2            = %9.4f \n',sige);
if results.rdraw == 0
fprintf(fid,'r-value            = %6d   \n',results.r);
else
fprintf(fid,'mean of rdraws     = %9.4f \n',mean(results.rdraw));
fprintf(fid,'gam(m,k) prior     = %6d,%6d \n',results.m,results.k);
end;    
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'# censored values  = %6d \n',results.nobsc);
fprintf(fid,'ndraws,nomit       = %6d,%6d \n',results.ndraw,results.nomit);
fprintf(fid,'accept rho rate    = %9.4f \n',results.acceptr);
fprintf(fid,'accept lam rate    = %9.4f \n',results.acceptl);
fprintf(fid,'time in secs       = %9.4f   \n',results.time);
fprintf(fid,'min and max rho    = %9.4f,%9.4f \n',results.rmin,results.rmax);
fprintf(fid,'min and max lambda = %9.4f,%9.4f \n',results.lmin,results.lmax);
fprintf(fid,'***************************************************************\n');

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
fprintf(fid,'      Posterior Estimates \n');

 if strcmp(results.pflag,'tstat')
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

% end of sact_g case

case {'sacp_g'} % <=================== Gibbs general spatial Probit model

nobs = results.nobs;
nvar = results.nvar;

% handling of vnames
Vname = 'Variable';
 for i=1:nvar
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
% add spatial rho parameter name
    Vname = strvcat(Vname,'lambda');

if (nflag == 1) % the user supplied variable names
Vname = 'Variable';
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_spat -- check vnames argument \n');
 fprintf(fid,'will use generic variable names \n');
 nflag = 0;
 else,
 for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
 end;
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    Vname = strvcat(Vname,'lambda');
 end; % end of if-else
end; % end of nflag issue

% find posterior means
tmp1 = mean(results.bdraw);
pout = mean(results.pdraw);
lout = mean(results.ldraw);
bout = [tmp1'
        pout
        lout];

y = results.y;
sige = mean(results.sdraw);
tmp1 = std(results.bdraw);
tmp2 = std(results.pdraw);
tmp3 = std(results.ldraw);
bstd = [tmp1'
        tmp2
        tmp3];

if strcmp(results.pflag,'tstat')
 tstat = bout./bstd;
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw results.ldraw];
 for i=1:results.nvar+2;
 if bout(i,1) > 0
 cnt = find(draws(:,i) > 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 else
 cnt = find(draws(:,i) < 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 end; % end of if - else
 end; % end of for loop
end; 

yhat = results.yhat;

fprintf(fid,'\n');
fprintf(fid,'Bayesian general spatial Probit model \n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'McFadden R^2    = %9.4f \n',results.r2mf);
fprintf(fid,'Estrella R^2    = %9.4f \n',results.rsqr);
fprintf(fid,'sigma^2         = %9.4f \n',sige);
if results.rdraw == 0
fprintf(fid,'r-value         = %6d   \n',results.r);
else
fprintf(fid,'mean of rdraws  = %9.4f \n',mean(results.rdraw));
fprintf(fid,'gam(m,k) prior  = %6d,%6d \n',results.m,results.k);
end; 
fprintf(fid,'Nobs, Nvars     = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'# 0, 1 y-values = %6d,%6d \n',results.zip,nobs-results.zip);
fprintf(fid,'ndraws,nomit    = %6d,%6d \n',results.ndraw,results.nomit);
fprintf(fid,'rho accept rate = %9.4f \n',results.acceptr);
fprintf(fid,'lam accept rate = %9.4f \n',results.acceptl);
fprintf(fid,'time in secs    = %9.4f   \n',results.time);
fprintf(fid,'min and max rho = %9.4f,%9.4f \n',results.rmin,results.rmax);
fprintf(fid,'***************************************************************\n');

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
fprintf(fid,'      Posterior Estimates \n');

 if strcmp(results.pflag,'tstat')
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
 
% <=================== end of sacp_g case

otherwise
error('results structure not known by prt_spat function');
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



