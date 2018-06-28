%addpath('C:\Users\BabaChattopadhyay\Desktop\Langmuir_journal\geom2d\geom2d')
function[volume]= pappus(out)
%theta=0:.00001:90;
%x=cosd(theta)';
%y=sind(theta)';
if (rem(length(out(:,1)),2)==0)
    
x=out(1:(length(out(:,1)))/2,1);
y=out(1:(length(out(:,1)))/2,2);
else
   x=out(1:((length(out(:,1)))-1)/2,1);
y=out(1:((length(out(:,1)))-1)/2,2); 
end
% x(length(x)+1)=0;
% y(length(y)+1)=0;
% X(:,1)=x;
% X(:,2)=y;
% PTS = centroid(X)
% A=polyarea(x,y);
% z=A*2*pi*PTS(1)
A=0;
area=0;
for i=1:length(x)-1
%     A=A+((x(i)+x(i+1))/2)*(polyarea([x(i),x(i),x(i+1),x(i+1)],[0,y(i),y(i+1),0]));
%     A=A+x(i)*(polyarea([x(i),x(i),x(i+1),x(i+1)],[0,y(i),y(i+1),0]));
   A=A+((x(i)+x(i+1))/2)*1/2*(y(i)+y(i+1))*(x(i)-x(i+1));
   area=area+1/2*(y(i)+y(i+1))*(x(i)-x(i+1));
end
cen=A/area;
vol=2*pi*cen*area;
volume=vol;
end