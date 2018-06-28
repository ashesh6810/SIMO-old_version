clear all;
close all;
clc;
N=30;


 theta=60;%% any angle of choice for a circle %%
R=2.5;
alpha=90-theta:2*theta/N:90+theta;
% 
pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];
% 

%
% pp(:,1)=[10;9;4;3;2;1;0];
% pp(:,2)=[0;1;7;5;4;1;0];



[out,normal_vec]=mvsplint(pp,N);


% theta=atand(out(1,3));
% if theta>0
%     theta=180-theta;
% end
% normal_vec=correction_fn(out,U,N,theta);

figure(2)
plot(out(:,1),out(:,2));
daspect([1 1 1])
hold on
% [normal_vec]=correction_fn(out,U,N,theta);
% disp(normal_vec);
 figure(2)
        quiver(out(:,1),out(:,2),normal_vec(:,1),normal_vec(:,2));
        daspect([1 1 1]);
        hold on;
%      original_area=polyarea(out(:,1),out(:,2));
%      [mn,min_ind]=min(out(2:N-1,4));
%      out(1,1) = out(1,1) + .1/norm(out(2,1:2)-out(1,1:2));   
%     out(N,1) = out(N,1) - .1/norm(out(N,1:2)-out(N-1,1:2));
%  
%     area_new=polyarea(out(:,1),out(:,2));
%     k=.001;
%     while(area_new-original_area>.0001)
%         out(min_ind,1:2)=out(min_ind,1:2)+k*normal_vec(min_ind,1:2);
%         
% 
% [out_updated,normal_vec]=mvsplint(out,N);
% out=out_updated;
% [minimum, min_ind] = min(out(2:N-1,4));
% figure(2)
% plot(out(:,1),out(:,2));
% hold on
% quiver(out(:,1),out(:,2),normal_vec(:,1),normal_vec(:,2));
% hold off
% 
%     end
    
       
    figure(2)
    plot(out(:,1),out(:,4),'r');
    hold on;
    figure(2)
    plot([0,10],[0,0],'b');
    
        
     