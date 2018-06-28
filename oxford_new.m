clear all;
close all;
p=1;
theta=60;
for i=0:.1:10
    x(p)=i;
    p=p+1;
end
for i=1:101
    y(i)=sqrt(((x(101)-x(1))/2*cscd(theta))^2-((x(i)-(x(101)+x(1))/2)^2))-((x(101)-x(1))/2)*cotd(theta);   
end
plot(x,y)
C=(x(101)+x(1))/2;
R=(x(101)-x(1))/2*cscd(theta);
A=R^2/2*((2*theta)*pi/180-sind(2*theta));

strip=.1; % can't change %
dA=.05;
dA_dT=.05;
M=1000;
theta1=60;
theta2=30;
tt=[];
pp(:,1)=x(1,:);
pp(:,2)=y(1,:);
out=mvsplint(pp,101);
count=1;
tt=[];
tt(count)=0;
area_arr=[];
ca=out(1,3)
for i=1:M
   
  
     GIF=floor((abs(out(101,1))-0)/(strip));
   remainder=rem(GIF,2);
   if remainder==0;
       h=theta2;
   else 
       h=theta1;
   end
    %% pinning %%
    if ca<h
          A=A+dA;
          w=1;
        zeta=ca;
        phi=.5;
        r=(x(101)-x(1))
        ashesh=1;
        while(ashesh==1||((A-A2)>.01))
           ashesh=2;
           
        zeta=zeta+phi;
        R=((x(101)-x(1))/2)*cscd(zeta)
        A2=R^2/2*((2*zeta)*pi/180-sind(2*zeta));
        for q=1:101
        y(q)=sqrt(((x(101)-x(1))/2*cscd(zeta))^2-((x(q)-(x(101)+x(1))/2)^2))-((x(101)-x(1))/2)*cotd(zeta);   
        end
        
    pp(:,1)=x(1,:);
    pp(:,2)=y(1,:);
    out1=mvsplint(pp,101);
     
        end
        ca=zeta
        pp(:,1)=out1(:,1);
        pp(:,2)=out1(:,2);
    end
    
        
    
    
    
        
        
     
    
  
   %% if advancing angle %% 
%     if ((out(1,3)-h)>0 && (out(1,3)-h)<.4)
if ca==h
        w=2;
        A=A+dA;
        R=sqrt(2*A/((2*h)*pi/180-sind(2*h)));
          
    r=R*sind(h);
    x(101)=r+C;
    x(1)=C-r;
    
    
    for q=1:101
        y(q)=sqrt(((x(101)-x(1))/2*cscd(h))^2-((x(q)-(x(101)+x(1))/2)^2))-((x(101)-x(1))/2)*cotd(h);   
    end
    pp(:,1)=x(1,:);
    pp(:,2)=y(1,:);
    
%     A=R^2/2*((2*h)*pi/180-sind(2*h));
end
   
     %% stretching %%
    if ca>h
        
    w=3;
    R=sqrt(2*A/((2*h)*pi/180-sind(2*h)));
          
    r=R*sind(h);
    x(101)=r+C;
    x(1)=C-r;
    
    
    for q=1:101
        y(q)=sqrt(((x(101)-x(1))/2*cscd(h))^2-((x(q)-(x(101)+x(1))/2)^2))-((x(101)-x(1))/2)*cotd(h);   
    end
    pp(:,1)=x(1,:);
    pp(:,2)=y(1,:);
    
    ca=h;

    end
    
   
     
    
    %% registering values %%
    out=mvsplint(pp,101);
    angle(count)=ca;
    radius(count)=R;
  
    energy(count)=R*2*ca*pi/180+cosd(ca)*(out(101,1)-out(1,1));
    area_arr(count)=A/(strip*10)^2;
    base_radius(count)=r;
    count=count+1;
    
    %% plotting %%
    figure(2)
    subplot(3,1,1)
    plot(out(:,1),out(:,2))
   axis([-10 20 0 20])
   axis equal
    hold off
    subplot(3,1,2)
    plot(area_arr,energy)
    axis([20 40 10 40])
    xlabel('S/a^2')
    ylabel('gibbs energy')
    hold on
    subplot(3,1,3)
    plot(area_arr,base_radius)
     axis([20 40 7 20])
     xlabel('S/a^2')
    ylabel('r/a^2')
    hold on
%     subplot(4,2,4)
%     plot(area_arr,angle)
%      axis([20 40 25 60])
%      xlabel('S/a^2')
%     ylabel('theta')
%     hold on
end