function F =myfun(x)
addpath('C:\Users\BabaChattopadhyay\Desktop\Langmuir_journal\inverse calculation');
X=xlsread('C:\Users\BabaChattopadhyay\Desktop\Langmuir_journal\inverse calculation\Book1.xlsx','Sheet1','A2:A29');
Y=xlsread('C:\Users\BabaChattopadhyay\Desktop\Langmuir_journal\inverse calculation\Book1.xlsx','Sheet1','B2:B29');
pp(:,1)=X;
pp(:,2)=Y;
N=60;
[out]=mvsplint(pp,N);


for k=1:N
F=x(1)-x(2)*out(k,2)-2*x(3)*out(k,4);
end
end















