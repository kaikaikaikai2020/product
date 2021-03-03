function [parameters, likelihood] = fattailed_garch2(data , p , q , breakpt , startingvals, options)
% PURPOSE:
%     FATTAILED_GARCH(P,Q) parameter estimation with different error distributions, the NOrmal, The T, 
%           and the Generalized Error Distribution
% 
% USAGE:
%     [parameters, likelihood, stderrors, robustSE, ht, scores] = fattailed_garch(data , p , q , errors, startingvals, options)
% 
% INPUTS:
%     data: A single column of zero mean random data, normal or not for quasi likelihood
% 
%     P: Non-negative, scalar integer representing a model order of the ARCH 
%       process
% 
%     Q: Positive, scalar integer representing a model order of the GARCH 
%       process: Q is the number of lags of the lagged conditional variances included
%       Can be empty([]) for ARCH process
% 
%     error:  The type of error being assumed, valid types are:
%            'NORMAL' - Gaussian Innovations
%            'STUDENTST' - T-distributed errors
%            'GED' - General Error Distribution
%  
%     startingvals: A (1+p+q) (plus 1 if STUDENTT OR GED is selected for the nu parameter) vector of starting vals.
%       If you do not provide, a naieve guess of 1/(2*max(p,q)+1) is used for the arch and garch parameters,
%       and omega is set to make the real unconditional variance equal
%       to the garch expectation of the expectation.
% 
%     options: default options are below.  You can provide an options vector.  See HELP OPTIMSET
% 
% OUTPUTS:
%     parameters : a [1+p+q X 1] column of parameters with omega, alpha1, alpha2, ..., alpha(p)
%                  beta1, beta2, ... beta(q)
% 
%     likelihood = the loglikelihood evaluated at he parameters
% 
%     robustSE = QuasiLikelihood std errors which are robust to some forms of misspecification(see White 94)
% 
%     stderrors = the inverse analytical hessian, not for quasi maximum liklihood
% 
%     ht = the estimated time varying VARIANCES
% 
%     scores = The numberical scores(# fo params by t) for M testing   
% 
% 
% COMMENTS:
%   GARCH(P,Q) the following(wrong) constratins are used(they are right for the (1,1) case or any Arch case
%     (1) Omega > 0
%     (2) Alpha(i) >= 0 for i = 1,2,...P
%     (3) Beta(i)  >= 0 for i = 1,2,...Q
%     (4) sum(Alpha(i) + Beta(j)) < 1 for i = 1,2,...P and j = 1,2,...Q
%     (5) nu>2 of Students T and nu>1 for GED
%
%   The time-conditional variance, H(t), of a GARCH(P,Q) process is modeled 
%   as follows:
%
%     H(t) = Omega + Alpha(1)*r_{t-1}^2 + Alpha(2)*r_{t-2}^2 +...+ Alpha(P)*r_{t-p}^2+...
%                    Beta(1)*H(t-1)+ Beta(2)*H(t-2)+...+ Beta(Q)*H(t-q)
%
%   Default Options
%   
%   options  =  optimset('fmincon');
%   options  =  optimset(options , 'TolFun'      , 1e-003);
%   options  =  optimset(options , 'Display'     , 'iter');
%   options  =  optimset(options , 'Diagnostics' , 'on');
%   options  =  optimset(options , 'LargeScale'  , 'off');
%   options  =  optimset(options , 'MaxFunEvals' , '400*numberOfVariables');
%
%
%  uses fFATTAILED_GARCHLIKELIHOOD and GARCHCORE.  You should MEX, mex 'path\garchcore.c', the MEX source 
%  The included MEX is for R12 Windows and was compiled with VC++6. It gives a 10-15 times speed improvement
%
%
% Author: Kevin Sheppard
% kksheppard@ucsd.edu
% Revision: 1    Date: 10/15/2000




t=size(data,1);
if nargin<6
    options=[];
end

errortype = 1;



if size(data,2) > 1
   error('Data series must be a column vector.')
elseif isempty(data)
   error('Data Series is Empty.')
end


if (length(q) > 1) | any(q < 0)
   error('Q must ba a single positive scalar or 0 for ARCH.')
end

if (length(p) > 1) | any(p <  0)
   error('P must be a single positive number.')
elseif isempty(p)
   error('P is empty.')
end

if isempty(q) | q==0;
   q=0;
   m=p;
else
   m  =  max(p,q);   
end


if nargin<=4 | isempty(startingvals)
   guess  = 1/(2*m+1);
   alpha  =  .15*ones(p,1)/p;
   beta   =  .75*ones(q,1)/q;
   omega  = (1-(sum(alpha)+sum(beta)))*cov(data(1:breakpt)) 
   omega2  =(1-(sum(alpha)+sum(beta)))*cov(data(breakpt+1:t))%set the uncond = to its expection
   nu=[];
  else
   omega=startingvals(1);
   omega2=startingvals(2);
   alpha=startingvals(3:p+2);
   beta=startingvals(p+3:p+q+2);
   nu=[];
   end
end


LB         =  [];     
UB         =  [];     
sumA =  [-eye(2+p+q); ...
      0  0 ones(1,p)  ones(1,q)];
sumB =  [zeros(2+p+q,1);...
      1];                          


if (nargin <= 5) | isempty(options)
   options  =  optimset('fmincon');
   options  =  optimset(options , 'TolFun'      , 1e-006);
   options  =  optimset(options , 'Display'     , 'iter');
   options  =  optimset(options , 'Diagnostics' , 'on');
   options  =  optimset(options , 'LargeScale'  , 'off');
   options  =  optimset(options , 'MaxFunEvals' , 400*(2+p+q));
end

sumB = sumB - [zeros(2+p+q,1); 1]*2*optimget(options, 'TolCon', 1e-6);

   LB = [];


   startingvals = [omega ;omega2; alpha ; beta];

% Estimate the parameters.
stdEstimate =  std(data(1:breakpt));  
stdEstimate2 =  std(data(breakpt+1:t));  
data    =  [stdEstimate(ones(m,1)) ; data];  
T=size(data,1);

[parameters, LLF, EXITFLAG, OUTPUT, LAMBDA, GRAD, HESSIAN] =  fmincon('fattailed_garchlikelihood2', startingvals ,sumA  , sumB ,[] , [] , LB , UB,[],options, data , p , q, errortype, stdEstimate, stdEstimate2, T, breakpt);


if EXITFLAG<=0
   EXITFLAG
   fprintf(1,'Not Sucessful! \n')
end

parameters(find(parameters    <  0)) = 0;          
parameters(find(parameters(1) <= 0)) = realmin;    
