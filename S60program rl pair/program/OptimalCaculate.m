function hOptimal=OptimalCaculate(sigma,mu,kappa,theta,eta,rho,lambda,T,t,x)
%% ����˵��
% �������Ž�h(t,x)
% �������lambdaΪ�������ϵ����TΪ��������Ϊ��λ
%% 
%lambda = 1/10;
alpha=kappa*(1-sqrt(1-lambda))/(2*eta^2)*(1+2*sqrt(1-lambda)/((1-sqrt(1-lambda))-(1+sqrt(1-lambda))*exp(2*kappa*(T-t)/sqrt(1-lambda))));
%beta=1/(2*eta^2*((1-sqrt(1-lambda))-(1+sqrt(1-lambda))*exp(2*kappa*(T-t)/sqrt(1-lambda))))*...
%    (lambda*sqrt(1-lambda)*(eta^2+2*rho*sigma*eta)*(1-exp(2*kappa*(T-t)/sqrt(1-lambda)))^2-...
%    lambda*(eta^2+2*rho*sigma*eta+2*kappa*theta)*(1-exp(2*kappa*(T-t)/sqrt(1-lambda))));
beta = kappa*theta*(1+sqrt(1-lambda))*(exp(2*kappa*(T-t)/sqrt(1-lambda))-1)/((1+(1-2/(1-sqrt(1-lambda)))*exp(2*kappa*(T-t)/sqrt(1-lambda))));
hOptimal=1/(1-lambda)*(beta+2*x*alpha-kappa*(x-theta)/eta^2+rho*sigma/eta+1/2);

end