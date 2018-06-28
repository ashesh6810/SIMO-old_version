clc
clear all;
close all;

N=100;
theta=10;
R=1;
alpha=90-theta:2*theta/N:90+theta;
pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];
count=1;
for i=1:10
    theta_counter(count)=theta
    alpha=90-theta:2*theta/N:90+theta;
    pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];
  vol1(count)= (pi/3)*(2-3*cosd(theta)+cosd(theta)^3)*R^3;
  vol2(count)=pappus(pp);
  theta=theta+5;
  count=count+1;
  
end
plot(theta_counter,vol1,'k');
hold on;
plot(theta_counter,vol2,'r*');
