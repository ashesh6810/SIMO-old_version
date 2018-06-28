close all;
clear all;
% num1= xlsread('C:\Users\BabaChattopadhyay\Desktop\rishi_codes\Book1.xlsx','Sheet1');
% num2= xlsread('C:\Users\BabaChattopadhyay\Desktop\rishi_codes\Book1.xlsx','Sheet2');
% figure(1)
% plot(num1,num2);
% hold on;

%% defining initial parameters %%
R=.1;
s=2;
div=10;
pitch=1/div;
N=25;
M=700;
Adot = 1;
theta=60;
dt = 0.004;
alpha=90-theta:2*theta/N:90+theta;

pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];

% out=mvsplint(pp,N);

theta_odd=60;
theta_even=30;
lca=theta;
RR=R;
 A=RR^2*(lca*pi/180-sind(2*lca)/2)/2;
 r=(pp(1,1)-pp(N,1))/2;


 %% start %%
 counter=1;
 tt=1;
 for count=1:M
    
     
   
    peri=2*lca*pi/180*RR;
%       for i=1:N-1
%     peri=peri+sqrt((out(i+1,1)-out(i,1))^2+(out(i+1,2)-out(i,2))^2);
%       end
    
    perimeter(counter)=peri;
    perturbation_factor=1000*N/peri;
% %     
%      perturbation_L = (R^2*(pi*lca/90 -sind(2*lca))/(perturbation_factor*sind(lca)));
 perturbation_L=.01;
     if mod(floor(div*pp(1,1)+1), s) == 0 
         theta_comp=theta_even;
     else
       theta_comp=theta_odd;
     end
if ((lca-theta_comp)==0)
    pp(1,1)=pp(1,1)+perturbation_L;
    pp(N,1)=pp(N,1)-perturbation_L;
    r=(pp(1,1)-pp(N,1))/2;
    RR=r*cscd(lca);
    A=RR^2*(2*lca*pi/180-sind(2*lca))/2;
    lca=theta_comp;
end
if(lca-theta_comp)>0
    pp(1,1)=pp(1,1)+perturbation_L;
    pp(N,1)=pp(N,1)-perturbation_L;
    r=(pp(1,1)-pp(N,1))/2;
    
 angle_value= non_liner_solver(pp,N,area(counter-1),lca);
 lca=angle_value;
 if (lca-theta_comp)<.1
     lca=theta_comp;
 end
 RR=r*cscd(lca);
end

if (lca-theta_comp)<0
    lca=lca+.1;
 r=(pp(1,1)-pp(N,1))/2;
RR=r*cscd(lca);
A=RR^2*(2*lca*pi/180-sind(2*lca))/2;
end

  %% calculating number of strips %%
    comp=(pp(1,1)-floor(pp(1,1)/pitch));
    if (comp==0||comp>pitch/2)
        num_strips=2*floor(pp(1,1)/pitch)-1;
    else  comp<=pitch/2;
        num_strips=2*floor(pp(1,1)/pitch)+1;
    end
%% how many %%
if mod(floor(div*pp(N,1)), s)==0
    num_strips_theta_oddl=floor(num_strips/2)+1;
   
else
    num_strips_theta_oddl=floor(num_strips/2);
end
num_strips_theta_evenl=num_strips- num_strips_theta_oddl;
%% calculating the last angle %%
if num_strips_theta_evenl>num_strips_theta_oddl
    last_angle=theta_odd;
else
    last_angle=theta_even;
end
%% energy calculation %%

Energy_circle(counter)=RR*2*lca*pi/(180)-num_strips_theta_evenl*cosd(theta_even)- num_strips_theta_oddl*cosd(theta_odd)-2*cosd(last_angle)*comp;
  angle(counter)=lca;
     area(counter)=A;
     radius(counter)=r;
     
counter=counter+1;
c(tt)=counter;
tt=tt+1;
figure(1)
% subplot(3,1,1)
% axis([28 62 0 1])
% xlabel('S/a^2');
% ylabel('theta')
% 
% plot(area,angle)
% subplot(3,1,2)
% plot(area,radius)
% 
% xlabel('S/a^2');
% ylabel('radius');
% subplot(3,1,3)


plot(area/pitch^2,angle,'r');
% axis([0 90 0 62]);
xlabel('S/a^2');
ylabel('energy');

 end
    
    