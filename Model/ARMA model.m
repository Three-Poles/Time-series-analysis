%ARMA模型：
clear all
%建模
Y=iddata(y);         %y为时间序列
%自相关性和互相关性
figure(1);
autocorr(y);
figure(2);
parcorr(y);
%na和nc是ARMA模型的阶次
na=2;               
nc=8;
m=armax(Y,[na nc])  %na和nc是模型阶次
%预测
S=5;                %S是预测步数
ff=[y;zeros(S,1)];
p=iddata(ff);
P=predict(m,p,S);   %通过模型进行预测
G=get(P);
PF=G.OutputData{1,1}(length(y)+1:length(y)+S,1)   %预测值