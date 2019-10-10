%ARIMA模型：
clear all
%建模
Y=diff(y);                %差分计算
Z=iddata(Y);
% 自相关性和互相关性
figure(1);
autocorr(Y);
figure(2);
parcorr(Y);
%p和q是ARMA模型的阶次
p=1;
q=8;
m=armax(Z,[p q])
%预测
S=5;                      %S是预测步数
ff=[Y;zeros(S,1)];
p=iddata(ff);
P=predict(m,p,S);         %通过模型进行预测
G=get(P);
PD=G.OutputData{1,1}(length(Y)+1:length(Y)+S,1);
D=[Y;PD];
X_D=cumsum([y(1);D]);
PF=X_D(length(y)+1:end)   %预测值