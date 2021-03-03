function prt_sp(results,vnames,fid)
% PURPOSE: Prints output using spatial regression results structures
%---------------------------------------------------
% USAGE: prt_sp(results,vnames,fid)
% Where: results = a structure returned by a spatial panel regression 
%        vnames  = an optional vector of variable names
%        fid     = optional file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%--------------------------------------------------- 
%  NOTES: e.g. vnames = strvcat('y','const','x1','x2');
%         e.g. fid = fopen('ols.out','wr');
%  use prt_sp(results,[],fid) to print to a file with no vnames               
% --------------------------------------------------
%  RETURNS: nothing, just prints the spatial panel regression results
% --------------------------------------------------

if ~isstruct(results)
 error('prt_sp requires structure argument');
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
 error('Wrong # of arguments to prt_sp');
end;

nvar = results.nvar;
nobs = results.nobs;

% handling of vnames
Vname = 'Variable';
 for i=1:nvar
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;

if (nflag == 1) % the user supplied variable names
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_sp -- check vnames argument \n');
 fprintf(fid,'will use generic variable names \n');
 nflag = 0;
 else,
Vname = 'Variable';
 for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
 end;
 end; % end of if-else
end; % end of nflag issue

fprintf(fid,'\n');

switch results.meth
case {'psem'}
fprintf(fid,'Pooled model with spatial error autocorrelation, no fixed effects \n');
case {'semsfe'}
fprintf(fid,'Pooled model with spatial error autocorrelation and spatial fixed effects \n');
case {'semtfe'}
fprintf(fid,'Pooled model with spatial error autocorrelation and time period fixed effects \n');
case {'semstfe'}
fprintf(fid,'Pooled model with spatial error autocorrelation, spatial and time period fixed effects \n');
case {'psar'}
fprintf(fid,'Pooled model with spatially lagged dependent variable, no fixed effects \n');
case {'sarsfe'}
fprintf(fid,'Pooled model with spatially lagged dependent variable and spatial fixed effects \n');
case {'sartfe'}
fprintf(fid,'Pooled model with spatially lagged dependent variable and time period fixed effects \n');
case {'sarstfe'}
fprintf(fid,'Pooled model with spatially lagged dependent variable, spatial and time period fixed effects \n');
otherwise
error('results structure not known by prt_sp function');
end % switch statement

switch results.meth

case {'psem','semsfe','semtfe','semstfe'} % <=================== spatial panel error models

% add spatial rho parameter name
    Vname = strvcat(Vname,'spat.aut.');
    
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared       = %9.4f   \n',results.rsqr);
fprintf(fid,'Rbar-squared    = %9.4f   \n',results.rbar);
fprintf(fid,'sigma^2         = %9.4f   \n',results.sige);
fprintf(fid,'log-likelihood  = %16.8g  \n',results.lik);
fprintf(fid,'Nobs,Nvar,TNvar = %6d,%6d,%6d  \n',results.nobs,results.nvar,results.tnvar);
fprintf(fid,'# iterations    = %6d     \n',results.iter);
fprintf(fid,'min and max rho = %9.4f,%9.4f \n',results.rmin,results.rmax);
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
fprintf(fid,'time for t-stats   = %9.4f \n',results.time3);
end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.liter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

fprintf(fid,'***************************************************************\n');

bout = [results.beta
        results.rho];
        
% <=================== end of sem case

case {'psar','sarsfe','sartfe','sarstfe'} % <=================== spatial panel autoregressive models

% add spatial rho parameter name
    Vname = strvcat(Vname,'W*dep.var.');
    
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared         = %9.4f   \n',results.rsqr);
fprintf(fid,'Rbar-squared      = %9.4f   \n',results.rbar);
fprintf(fid,'sigma^2            = %9.4f \n',results.sige);
fprintf(fid,'Nobs,Nvar,TNvar    = %6d,%6d,%6d  \n',results.nobs,results.nvar,results.tnvar);
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
fprintf(fid,'iter  for MC appr  = %6d  \n',results.liter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

fprintf(fid,'***************************************************************\n');

bout = [results.beta
        results.rho];
        
% <=================== end of sar case

otherwise
error('results structure not known by prt_sp function');
end;

% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities, function of LeSage
tmp = [bout results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in); %function of LeSage