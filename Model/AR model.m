%AR模型:
clear all
%建模
Y=iddata(y);         %y为时间序列
figure(1);
parcorr(y);
n=2;                 %n是AR模型的阶次
m=ar(Y,n,'ls')       %n是AR的阶次，’ls’是参数估计方法，为最小二乘方法
%预测
S=5;                 %S是预测步数
ff=[y;zeros(S,1)];    
p=iddata(ff);
P=predict(m,p,S);    %通过模型进行预测
G=get(P);
PF=G.OutputData{1,1}(length(y)+1:length(y)+S,1)     %预测值