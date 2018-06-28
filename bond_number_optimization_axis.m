clear all;
close all;
N=25;
theta=60;
R=1*2/1.732;
Ro=1;
alpha=90-theta:2*theta/N:90+theta;
pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];
Bo=100;

M=10000;
[out,normal_vec]=mvsplint(pp,N);
outvol=mvsplint(out,1000);
vol=pappus(outvol);
%area=polyarea(out(:,1),out(:,2));

dV=.02*vol/100;

%%
% filename= [ 'Bond_number=',num2str(.13*R^2),'.avi'];     
% vidObj = VideoWriter(filename);
% open(vidObj);
%%
for count=1:M
    ashesh=0;
    %%
    peri=0;
    % Calculation of the perimeter
    for i=1:N-1
    peri=peri+sqrt((out(i+1,1)-out(i,1))^2+(out(i+1,2)-out(i,2))^2);
    end
%    perimeter(counter)=peri;
    % Deciding the perturbation factor for moving the end points
    perturbation_factor = 1000*N/peri;
    %%
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
    
     perturbation_R =( R^2*(pi*rca/90 -sind(2*rca))/(perturbation_factor*sind(rca)));
       perturbation_L = (R^2*(pi*lca/90 -sind(2*lca))/(perturbation_factor*sind(lca)));
    
    
    
    %%
    if lca>theta
        ashesh=1;
        outnew=stretch_droplet_Bo_axis(out, N, .0001, vol, 'both',normal_vec,Bo/(Ro)^2);
    else
          ashesh=2;
          outnew=pin_droplet_reced_Bo_axis(out,N,vol-dV,normal_vec,Bo/(Ro)^2);
        
    [outnew]=pin_droplet_Bo_axis(outnew,N,vol,normal_vec,Bo/(Ro)^2);
    
     %area=polyarea(outnew(:,1),outnew(:,2));
%     %dA_comp=area_new-area;

   end
    [out,normal_vec]=mvsplint(outnew,N);
    comp=(out(2:N-1,4)+out(2:N-1,8))-Bo/R^2*out(2:N-1,2);
%     outvolume=mvsplint(out,100);
%     vol=pappus(outvolume);
%     dV=vol-.99998*vol;
h=figure(1)
subplot(2,1,1)
plot(pp(:,1),pp(:,2),'k-');
daspect([1 1 1]);
hold on

    plot(out(:,1),out(:,2));
    daspect([1 1 1]);
    %axis([-5 5 0 10])
 hold off
 subplot(2,1,2)
 plot(comp)

%       F = getframe(h);
%    writeVideo(vidObj,F);
end