
function[radius,angle,area,Energy_circle]=kusumatmaja
%% defining initial parameters %%
R=1*2/1.732;
% s=2;
% div=10;
% pitch=1/div;
length_strip=1;

N=25;
MM=1700;
Adot = 1;
theta=60;
dt = 0.004;
alpha=90-theta:2*theta/N:90+theta;
area=[];
radius=[];
angle=[];

pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];

% out=mvsplint(pp,N);

theta_white=30;
theta_black=60;
lca=theta;
RR=R;
 A=RR^2*(lca*pi/180-sind(2*lca)/2)/2;
 r=(pp(1,1)-pp(N,1))/2;


 %% start %%
 counter=1;
 tt=1;
 for count=1:MM
    
     
   
    peri=2*lca*pi/180*RR;

    
    perimeter(counter)=peri;
    perturbation_factor=1000*N/peri;
% %     
%      perturbation_L = (R^2*(pi*lca/90 -sind(2*lca))/(perturbation_factor*sind(lca)));
 perturbation_L=.01;

%% deciding which strip %%
 I=(floor((pp(1,1)-(length_strip)/2)));
 uu=floor((pp(1,1)-(length_strip)/2)/2);
 if rem(I,2)==0
     theta_comp=theta_black;
 else
     theta_comp=theta_white;
 end
 
 %% calculating number of strips %%

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


%% energy calculation %%

Energy_circle(counter)=RR*2*lca*pi/(180)-num_black*cosd(theta_black)- num_white*cosd(theta_white)-cosd(theta_comp)*last_fraction;
  angle(counter)=lca;
  f(counter)=2*(lca/sin(lca)-cos(lca));
     area(counter)=A;
     radius(counter)=r;
     
counter=counter+1;
c(tt)=counter;
tt=tt+1;


 end
end
    
    