% PURPOSE: A Monte Carlo comparison of GMM sem and ML sem
% estimation of the spatial error model
% in a Monte Carlo experiment that compares maximum likelihood and GMM estimation  
%---------------------------------------------------
% USAGE: sar_mcarlo
%---------------------------------------------------

clear all;
n = 500; % sample size (nobs)

xc = randn(n,1); % random location coordinates
yc = randn(n,1);
% W-matrix using xy2cont
% create sparse standardized 1st-order spatial weight matrix
[j1 W j2] = xy2cont(xc,yc);

% compute eigenvalue range for rho
opt.tol = 1e-3; opt.disp = 0;
lambda = eigs(sparse(W),speye(n),1,'SR',opt);  
rmin = 1/lambda;   
rmax = 1;

% these are the true values for sige and beta 
sige = 0.5;
k = 3;
x = randn(n,k);
beta = zeros(k,1);
beta(1,1) = 1.0;
beta(2,1) = 1.0;
beta(3,1) = 1.0;

% do a grid of rho values
% [-0.8, -0.4, 0 0.4, 0.8]
rgrid = -0.8:0.4:0.8;

niter = 500; % # of experiments to carry out

nrho = length(rgrid);

% ========== begin the fun ==========================
for i=1:nrho; % loop over rho values
    rho = rgrid(i);
    
bout = zeros(niter,2*k); % storage for results
pout = zeros(niter,2);
sout = zeros(niter,2);


time1 = 0;
time2 = 0;

% ==========> begin looping
% we loop for niter times generating a new y-vector every time through
% keeping the x-matrix fixed
for iter=1:niter;

% generate a y, with fixed X-matrix and new epsilon vector
y = (speye(n) - rho*W)\x*beta + (speye(n) - rho*W)\(randn(n,1)*sqrt(sige));;

if iter == 1 % we speed things up by only computing the log-determinant
             % that stays the same once on the 1st iteration
info.lflag = 0;
info.rmin = rmin;
info.rmax = rmax;
results0 = sar(y,x,W,info);
lndetv = results0.lndet;
results1 = sar_gmm(y,x,W);
else
options.lndet = lndetv;
options.rmin = rmin;
options.rmax = rmax;
results0 = sar(y,x,W,options);
results1 = sar_gmm(y,x,W);
end;

% extract and save the results for a pretty print-out
bout(iter,1:k) = results0.beta';
sout(iter,1) = results0.sige;
pout(iter,1) = results0.rho;
time1 = time1 + results0.time;

bout(iter,k+1:2*k) = results1.beta';
sout(iter,2) = results1.sige;
pout(iter,2) = results1.rho;
time2 = time2 + results1.time;

end;
% ==========> end of looping

% compute things to make a results table

bounds1 = zeros(3,2);
bounds2 = zeros(3,2);

% find 95 percentile points
for i=1:3;
bounds1(i,:) = hpdi(bout(:,i),0.95);
bounds2(i,:) = hpdi(bout(:,k+i),0.95);
end;


% row and column labels for a pretty print out of the results

in.rnames = strvcat('Estimates','ML b1','ML b2','ML b3','GM b1','GM b2','GM b3', ...
    'ML rho','GMM rho','ML sige','GMM sige');
in.cnames = strvcat('min','0.05','means','median','Truth','0.95','max');
in.width = 1000;

bprint = [min(bout)
          bounds1(:,1)' bounds2(:,1)' 
          mean(bout)
          median(bout)
          [beta' beta']
          bounds1(:,2)' bounds2(:,2)' 
          max(bout)];


bbounds1 = zeros(2,2);
bbounds2 = zeros(2,2);

for i=1:2
bbounds1(i,:) = hpdi(pout(:,i),0.95);
bbounds2(i,:) = hpdi(sout(:,i),0.95);
end;

pprint1 = [min(pout)  
          bbounds1(:,1)'
          mean(pout)  
          median(pout)
          [rho rho]
          bbounds1(:,2)'
          max(pout)];
pprint2 = [min(sout)  
          bbounds2(:,1)'
          mean(sout)  
          median(sout)
          [sige sige]
          bbounds2(:,2)'
          max(sout)];
                  

fprintf(1,'============================================================ \n');
fprintf(1,'results from %5d Monte Carlo simulations rho = %5.2f, nobs = %10d \n',niter,rho,n);

out = [bprint'
       pprint1'
       pprint2'];
   
mprint(out,in);

in3.cnames = strvcat('ML TIME','GMM TIME');
mprint([time1 time2],in3);

end; % end of loop over different rho values
