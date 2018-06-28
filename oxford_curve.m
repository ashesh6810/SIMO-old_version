clear all;
close all;
p=1;
theta=60;
N=100;
alpha=90-theta:2*theta/N:90+theta;
R=1;
pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];


plot(pp(:,1),pp(:,2))
daspect([1 1 1])

out=mvsplint(pp,N);
strip=1;
theta1=60;
theta2=30;
dA_dT=1;
area=R^2/2*(180-2*theta)*pi/180;
perimeter=R*pi/180*(180-theta);
E=perimeter+cosd(theta)*(out(N,1)-out(1,1));
epsilon=.5;


M=1000;

for i=1:M
    area=area+dA_dT;
    remainder=rem(out(1,1),2);
    if remainder==0
        h=theta1;
    else h=theta2;
    end
    if out(1,3)>h
        out(1,1)=out(1,1)+epsilon;
        out(N,1)=out(N,1)-epsilon;
       R=sqrt(2*area/(pi/180*(180-h)));
        alpha=90-theta:2*theta/N:90+theta;
        pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];
        out=mvsplint(pp,N);
        aaaa=1;
    else
        R=sqrt(2*area/(pi/180*(180-h)));
        alpha=90-theta:2*theta/N:90+theta;
        pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];
        out=mvsplint(pp,N);
        aaaa=2;
    
    end
    area=R^2/2*(180-2*theta)*pi/180;
perimeter=R*pi/180*(180-theta);
    plot(out(:,1),out(:,2));
    hold off
    
end

    
    
        
        
        
    