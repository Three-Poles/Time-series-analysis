%MA模型：
clear all
%建模
Y=iddata(y);        %y为时间序列
figure(1);
autocorr(y);
n=8;                %n是MA模型的阶次
m=armax(Y,'nc',n)   %n是MA的阶次
%预测
S=5;                %S是预测步数
ff=[y;zeros(S,1)];
p=iddata(ff);
P=predict(m,p,S);   %通过模型进行预测
G=get(P);
PF=G.OutputData{1,1}(length(y)+1:length(y)+S,1)    %预测值