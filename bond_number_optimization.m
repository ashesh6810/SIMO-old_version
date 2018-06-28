clear all;
close all;
N=25;
theta=60;
R=5;
alpha=90-theta:2*theta/N:90+theta;
pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];

Bo=5;
M=10000;
[out,normal_vec]=mvsplint(pp,N);
area=polyarea(out(:,1),out(:,2));

dA=area-.9998*area;
%%
filename= [ 'Bond_number=',num2str(Bo/R^2),'.avi'];     
vidObj = VideoWriter(filename);
open(vidObj);
%%
for count=1:M
    %% contact angle
    if out(N,3)<=0
    lca = 180+(out(N,3));
    else 
        lca=out(N,3);
    end
    
    if out(1,3)>=0;
    rca = 180-out(1,3);
    else
        rca =abs(out(1,3) );                      %abs to eliminate negative - will not work for theta>90
    end
    %%
    if lca>theta
        outnew=stretch_droplet_Bo(out, N, .001, area, 'both',normal_vec,Bo);
    else
    [outnew]=pin_droplet_Bo(out,N,area+dA,normal_vec,Bo/(2*R)^2);
    
     %area=polyarea(outnew(:,1),outnew(:,2));
%     %dA_comp=area_new-area;
  outnew=pin_droplet_reced_Bo(outnew,N,area,normal_vec,Bo/R^2);
    end
    [out,normal_vec]=mvsplint(outnew,N);
h=figure(1)
    plot(outnew(:,1),outnew(:,2));
    daspect([1 1 1]);
    %axis([-5 5 0 10])
%   hold off
      F = getframe(h);
   writeVideo(vidObj,F);
end