function [sigma,mu,kappa,theta,eta,rho]=ParaEstimate(St,Xt,N,dt,minKappa)
%% 函数说明
% 模型涉及的参数估计
% 输入变量包括St=ln(Bt),Xt=ln(At)-ln(Bt),N为学习时间,dt为观察时间间隔,其中St和Xt均为N+1向量
%%
mIn=(St(N+1)-St(1))/N;
SIn=sqrt((sum((St(2:N+1)-St(1:N)).^2)-2*mIn*(St(N+1)-St(1))+N*mIn^2)/N);
pIn=1/(N*sum(Xt(1:N).^2)-sum(Xt(1:N))^2)*(N*sum(Xt(2:N+1).*Xt(1:N))-(Xt(N+1)-Xt(1))*sum(Xt(1:N))-sum(Xt(1:N))^2);
qIn=(Xt(N+1)-Xt(1)+sum(Xt(1:N))-pIn*sum(Xt(1:N)))/N;
VIn=sqrt(abs(1/N*(Xt(N+1)^2-Xt(1)^2+(1+pIn^2)*sum(Xt(1:N).^2)-2*pIn*sum(Xt(1:N).*Xt(2:N+1))-N*qIn)));
CIn=1/(N*VIn*SIn)*(sum(Xt(2:N+1).*(St(2:N+1)-St(1:N)))-pIn*sum(Xt(1:N).*(St(2:N+1)-St(1:N)))-mIn*(Xt(N+1)-Xt(1))-mIn*(1-pIn)*sum(Xt(1:N)));

sigma=sqrt(SIn^2/dt);
mu=mIn/dt+1/2*sigma^2;
kappa=-log(pIn)/dt;
theta=qIn/(1-pIn);
eta=sqrt(2*kappa*VIn^2/(1-pIn^2));
rho=kappa*CIn*VIn*SIn/(eta*sigma*(1-pIn));

% 参数优化
kappa=max(kappa,minKappa);
end