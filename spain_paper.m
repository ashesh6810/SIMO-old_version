 close all;
clear all;
% num1= xlsread('C:\Users\BabaChattopadhyay\Desktop\rishi_codes\Book1.xlsx','Sheet1');
% num2= xlsread('C:\Users\BabaChattopadhyay\Desktop\rishi_codes\Book1.xlsx','Sheet2');
% figure(1)
% plot(num1,num2);
% hold on;

%% defining initial parameters %%
%R=1*2/1.732;
R=.2;
% s=2;
% div=10;
% pitch=1/div;
length_strip=.1;

N=25;
M=10000;
Adot = 1;
theta=57.5;
dt = 0.004;
alpha=90-theta:2*theta/N:90+theta;

pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];

% out=mvsplint(pp,N);

theta_white=57.5;
theta_black=128;
lca=theta;
RR=R;
 %A=RR^2*(lca*pi/180-sind(2*lca)/2)/2;
 A=(RR^3/3)*pi*(2-3*cosd(lca)+cosd(lca)^3);
 r=(pp(1,1)-pp(N,1))/2;
%%
 filename= ['kusumatmaja','.avi'];     
 vidObj = VideoWriter(filename);
 open(vidObj);

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
%      if mod(floor(div*pp(1,1)+1), s) == 0 
%          theta_comp=theta_even;
%      else
%        theta_comp=theta_odd;
%      end
%% deciding which strip %%
 I=(floor((pp(1,1)-(length_strip)/2)));
 uu=floor((pp(1,1)-(length_strip)/2)/2);
 if rem(I,2)==0
     theta_comp=theta_black;
 else
     theta_comp=theta_white;
 end
 
 %% calculating number of strips %%
%  II=floor((pp(1,1)-.5)/2);
%  if rem(II,2)==0
%  num_white=2*floor(((pp(1,1)-.5)/2))+1;
% num_black=2*floor(((pp(1,1)-.5)/2));
%  else
%      num_white=2*floor(((pp(1,1)-.5)/2))+1;
% num_black=2*floor(((pp(1,1)-.5)/2))+2;
%  end
if rem(I,2)==0
    
        num_white=2*uu+1;
        num_black=2*uu;
else
          num_white=2*uu+1;
    num_black=2*uu+2;
end
        

%% calculating last fraction %%
last_fraction =2*((pp(1,1)-.5)-2*floor((pp(1,1)-.5)/2));

if ((lca-theta_comp)==0)
    pp(1,1)=pp(1,1)+perturbation_L;
    pp(N,1)=pp(N,1)-perturbation_L;
    r=(pp(1,1)-pp(N,1))/2;
    RR=r*cscd(lca);
   % A=RR^2*(2*lca*pi/180-sind(2*lca))/2;
   A=(RR^3/3)*pi*(2-3*cosd(lca)+cosd(lca)^3);
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
%A=RR^2*(2*lca*pi/180-sind(2*lca))/2;
A=(RR^3/3)*pi*(2-3*cosd(lca)+cosd(lca)^3);
end

  %% calculating number of strips %%
%     comp=(pp(1,1)-floor(pp(1,1)/pitch));
%     if (comp==0||comp>pitch/2)
%         num_strips=2*floor(pp(1,1)/pitch)-1;
%     else  comp<=pitch/2;
%         num_strips=2*floor(pp(1,1)/pitch)+1;
%     end


%% how many %%
% if mod(floor(div*pp(N,1)), s)==0
%     num_strips_theta_oddl=floor(num_strips/2)+1;
%    
% else
%     num_strips_theta_oddl=floor(num_strips/2);
% end
% num_strips_theta_evenl=num_strips- num_strips_theta_oddl;
%% calculating the last angle %%
% if num_strips_theta_evenl>num_strips_theta_oddl
%     last_angle=theta_odd;
% else
%     last_angle=theta_even;
% end
%% energy calculation %%

Energy_circle(counter)=RR*2*lca*pi/(180)-num_black*cosd(theta_black)- num_white*cosd(theta_white)-cosd(theta_comp)*last_fraction;
  angle(counter)=lca;
  f(counter)=2*(lca/sin(lca)-cos(lca));
     area(counter)=A;
     radius(counter)=r;
     
counter=counter+1;
c(tt)=counter;
tt=tt+1;

h=figure(1)

% axis([28 62 0 1])
plot((area.^(1/3)),angle)
xlabel('S/a^2');
ylabel('theta')

   F = getframe(h);
    writeVideo(vidObj,F);
 end
    close(vidObj);
    