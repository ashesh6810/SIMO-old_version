clear all;
clc;
close all;
addpath('C:\Users\BabaChattopadhyay\Desktop\Langmuir_journal\inverse calculation');
X=xlsread('C:\Users\BabaChattopadhyay\Desktop\Langmuir_journal\inverse calculation\Book2.xlsx','Sheet1','A1:A23');
Y=xlsread('C:\Users\BabaChattopadhyay\Desktop\Langmuir_journal\inverse calculation\Book2.xlsx','Sheet1','B1:B23');
pp(:,1)=X;
pp(:,2)=Y;
N=60;
[out]=mvsplint(pp,N);
for i=1:2:N
    out1(i,:)=out(i,:);
end    
%x0=[1,1,1];
%[x,resnorm,res,eflag,output1]=lsqnonlin(@myfun,x0)
for i=1:N
   C(i,1)=-2*(out(i,4)+out(i,8));
   C(i,2)=-1;
end
for i=1:N
d(1,i)=out(i,2);
end
x=lsqlin(C,d);
V=pappus(out);
a=(3*V/(4*pi))^(1/3);
Bo=(1/x(1))*(a)^2