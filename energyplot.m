clear all;
close all;

% plt = 1 if you want to see the plots as the code is ruuning
plt = 1;
% dA/dt
Adot = 1;
dt = 0.004;
% Intial contact angle of the droplet 
theta = 60; 
% Intial major radius of the droplet
R = 1;
% Pitch of the hetrogenity is 1/div
div = 10;
pitch=1/div;
s = 2;
% Number of points on the curve
N = 25;
% Number of loops to run
M=1700;

% IMM = imread('Picture.jpg');


%% Initial Calculations

% Define the initial shape
alpha=90-theta:2*theta/N:90+theta;

pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];

% The depinning contact angles of the right side
theta_oddr = 60;
theta_evenr = 30;

% The depinning contact angles of the left side
theta_oddl = 60;
theta_evenl = 30;

% Creating the spline
[outnew,normal_vec]=mvsplint(pp,N);
out = outnew;

% Calculating original area
area_orig = polyarea(pp(:,1),pp(:,2));
area_first = area_orig;

    
% Creating an array to define CAs
lca_arr=zeros(1,M);
rca_arr=zeros(1,M);
% area_arr=zeros(1,M);
% lca_circle=zeros(1,M);
ca_fitted=zeros(1,M);
t = zeros(1,M);
tt = zeros(1,M);
countt=1;

    % Calculating the two CAs
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
lca_circle=theta;
RR=R;
A=RR^2*(lca*pi/180-sind(2*lca)/2)/2;
r=(pp(1,1)-pp(N,1))/2;
% Clock starter
tic;
% filename= [num2str(round(rand()*1000)), 'Video_theta_1=', num2str(theta_oddl),'_theta_2=', num2str(theta_evenl),'.avi'];     
% vidObj = VideoWriter(filename);
% open(vidObj);
counter=1;
for count = 1:M   
    
    % Update area
    area_orig = polyarea(out(:,1),out(:,2)); 
    % Intialization of the perimeter
    peri=0;
    % Calculation of the perimeter
    for i=1:N-1
    peri=peri+sqrt((out(i+1,1)-out(i,1))^2+(out(i+1,2)-out(i,2))^2);
    end
    perimeter(counter)=peri;
    % Deciding the perturbation factor for moving the end points
    perturbation_factor = 1000*N/peri;
    

    
    % rightmost point is on even strip//1 added for symmetry//always negative
    if mod(floor(div*out(1,1)+1), s) == 0 
        delta_R = rca - theta_evenr;
        
    else
        delta_R = rca - theta_oddr;
        
    end
    %% calculating number of strips %%
    comp=(out(1,1)-floor(out(1,1)/pitch));
    if (comp==0||comp>pitch/2)
        num_strips=2*floor(out(1,1)/pitch)-1;
    else  comp<=pitch/2;
        num_strips=2*floor(out(1,1)/pitch)+1;
    end
    %%
    
    perturbation_R =( R^2*(pi*rca/90 -sind(2*rca))/(perturbation_factor*sind(rca)));
    
    % leftmost point is on even strip
    if mod(floor(div*out(N,1)), s)==0
        delta_L = lca - theta_evenl;
        
    else
        delta_L = lca - theta_oddl;
       
    end
    perturbation_L = (R^2*(pi*lca/90 -sind(2*lca))/(perturbation_factor*sind(lca)));
%% how many %%
if mod(floor(div*out(N,1)), s)==0
    num_strips_theta_oddl=floor(num_strips/2)+1;
   
else
    num_strips_theta_oddl=floor(num_strips/2);
end
num_strips_theta_evenl=num_strips- num_strips_theta_oddl;
%% last angle %%
if num_strips_theta_evenl>num_strips_theta_oddl
    last_angle=theta_oddl;
else
    last_angle=theta_evenl;
end
%% energy %%
Energy(counter)=perimeter(counter)-num_strips_theta_evenl*cosd(theta_evenl)- num_strips_theta_oddl*cosd(theta_oddl)-2*cosd(last_angle)*comp;
    % decide which side to move
    if delta_R>=0 && delta_L>=0
        outnew = stretch_droplet(out, N, (perturbation_L+perturbation_R)/2, area_orig, 'both',normal_vec);
%         angle(counter)=non_liner_solver(outnew,N,area_orig,lca)

    else
        outnew = pin_droplet(out, N, area_orig + Adot*dt,normal_vec);
%         angle(counter)=non_liner_solver1(outnew,N,area_orig+Adot*dt,lca)
    end
    [out normal_vec] = mvsplint(outnew,N);
    rr(counter)=(out(1,1)-out(N,1))/2;
   %% circle validation %%
    if mod(floor(div*out(1,1)+1), s) == 0 
         theta_comp=theta_evenl;
     else
       theta_comp=theta_oddl;
    end
    if ((lca_circle-theta_comp)==0)
    pp(1,1)=pp(1,1)+perturbation_L;
    pp(N,1)=pp(N,1)-perturbation_L;
    r=(pp(1,1)-pp(N,1))/2;
    RR=r*cscd(lca_circle);
    A=RR^2*(2*lca_circle*pi/180-sind(2*lca_circle))/2;
    lca_circle=theta_comp;
    end
if(lca_circle-theta_comp)>0
    pp(1,1)=pp(1,1)+perturbation_L;
    pp(N,1)=pp(N,1)-perturbation_L;
    r=(pp(1,1)-pp(N,1))/2;
    
 angle_value= non_liner_solver(pp,N,area(counter-1),lca);
 lca_circle=angle_value;
 if (lca_circle-theta_comp)<.1
     lca_circle=theta_comp;
 end
 RR=r*cscd(lca_circle);
end
if (lca_circle-theta_comp)<0
    lca_circle=lca_circle+.1;
 r=(pp(1,1)-pp(N,1))/2;
RR=r*cscd(lca);
A=RR^2*(2*lca_circle*pi/180-sind(2*lca))/2;
end
 angle(counter)=lca_circle;
     area(counter)=A;
          

     %%
% Energy_circle(counter)=r(counter)*2*lca*pi/(180*sind(angle(counter)))-num_strips_theta_evenl*cosd(theta_evenl)- num_strips_theta_oddl*cosd(theta_oddl)-2*cosd(angle(counter))*comp;
    % calculate the current areas
    area_now  = polyarea(out(:,1),out(:,2));
     qq=1;
    y_circle=[];
    pp_circle_fit=[];
    
    for jj=1:N
     y_circle(qq)=sqrt(((out(N,1)-out(1,1))/2*cscd(angle(counter)))^2-((out(jj,1)-(out(N,1)+out(1,1))/2)^2))-((out(N,1)-out(1,1))/2)*cotd(angle(counter));
     qq=qq+1;
    end
    
   
%     circleout=mvsplint(pp_circle_fit,N);
     
    
      %area of circle throught the end points %
         area_arr(counter)=area_now;

    % fitting circle through end points %
    
   counter=counter+1;
    % calculate time spent
    t(count) = (area_now-area_first)/Adot;
    
    %% test plot %%
%         figure(11)
%         plot(pp_circle_fit(:,1),pp_circle_fit(:,2))
% if counter>100
% figure(22)
        plot(area,angle)
% %         
% %      axis([.6 2 38 70])
      hold on

        %%

    if count>1
    if t(count)>t(count-1)
        
        % fitting a circle to the shape
        [xc,yc, Re] = circfit(out(:,1),out(:,2));
        xfit = Re*cosd(0:.1:180)+xc; yfit = Re*sind(0:.1:180)+yc;
        ydash = sqrt(Re^2-yc^2)/yc;
        
      
         
        %ca_fitted(countt) = abs(atand(ydash));
        if out(N,3)<=0
            lca_arr(countt)= 180+(out(N,3))+4;
        else
            lca_arr(countt)=out(N,3)+4;
        end
        if atand((ydash))>=0
            ca_fitted(countt)= 180-(atand((ydash)));
        else
            ca_fitted(countt)= abs(atand((ydash)));
        end
        if out(1,3)>=0;
            rca_arr(countt) = 180-out(1,3)+4;
         else
            rca_arr(countt) =abs(out(1,3))+4;                        %abs to eliminate negative - will not work for theta>90
        end
         lca = lca_arr(countt);
         rca = rca_arr(countt);
        tt(countt) = t(count);
        
        %% fitting circle %%
       
        
        %%
        % if we want to plot
        if plt==1 
            h = figure(1);
            set(h, 'Position', [50 50 1024 640])
            subplot(10,3,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21]);
            plot(pp(:,1), pp(:,2), 'k', 'LineWidth',1);
            hold on;
            plot(out(:,1), out(:,2), 'LineWidth',1.5);
            hold on;
            axis([-4*R 4*R -R/4 2.75*R]);
            for q = 0:s:60
            plot([(q)/div (s+q-1)/div], [0 0], '-c', 'LineWidth',2);
            plot([(q+s-1)/div (q+s)/div], [0 0], '-r', 'LineWidth',2);
            plot([-(s+q-1)/div -(q)/div], [0 0], '-c', 'LineWidth',2);
            plot([-(q+s)/div -(s+q-1)/div], [0 0], '-r', 'LineWidth',2);
            end
            if tt(countt)>=4 && tt(countt)<6
            quiver(out(:,1), out(:,2),normal_vec(:,1).*abs(out(:,4)), normal_vec(:,2).*abs(out(:,4)), '-m', 'LineWidth',1.5);
            end
            if tt(countt)>=2 && tt(countt)<4
                if rca_arr(countt)<=90
                quiver(out(1,1), out(1,2),-abs(normal_vec(1,2)), abs(normal_vec(1,1)), '-m', 'LineWidth',1.5);
                quiver(out(N,1), out(N,2),abs(normal_vec(N,2)), abs(normal_vec(N,1)), '-m', 'LineWidth',1.5);
                end
                if rca_arr(countt)>90
                quiver(out(1,1), out(1,2),abs(normal_vec(1,2)), abs(normal_vec(1,1)), '-m', 'LineWidth',1.5);
                quiver(out(N,1), out(N,2),-abs(normal_vec(N,2)), abs(normal_vec(N,1)), '-m', 'LineWidth',1.5);
                end
            end

            
           % axis off;
            daspect([1 1 1]);
            legend('Initial Shape','Current Shape', ['\theta_1 = ', num2str(theta_oddl), '^o'], ['\theta_2 = ', num2str(theta_evenl), '^o']);
            set(gcf,'Color',[1,1,1]);
            hold off;
            subplot(10,3,[23 24 26 27 29 30]);
            hold on;
            plot(tt(1:countt), lca_arr(1:countt),'-r', 'LineWidth',2);
            plot(tt(1:countt), rca_arr(1:countt),'-g', 'LineWidth',2);
            plot(tt(1:countt), ca_fitted(1:countt),'--k', 'LineWidth',2);
            axis([0 8 theta_evenl-5 theta_oddl+5]);
            legend('Actual Left','Actual Right', 'Apparent');
            set(gcf,'Color',[1,1,1]);
            xlabel('Time', 'FontSize', 12,'Color', [0 0 0]);
            ylabel('Contact Angle ( ^o )', 'FontSize', 12,'Color', [0 0 0]);
            grid on;
            hold off;
            
            subplot(10,3,[22 25 28]);
%             imshow(IMM);
            hold off;
        end
    countt = countt+1;
    F = getframe(h);
%     writeVideo(vidObj,F);
    end
    end
%      if countt==104
%          Ans1(1,1:25) = out(1:25,1);
%          Ans1(2,1:25) = out(1:25,2);
%      end
%      if countt==157
%          Ans2(1,1:25) = out(1:25,1);
%          Ans2(2,1:25) = out(1:25,2);
%      end
%      if countt==380
%          Ans3(1,1:25) = out(1:25,1);
%          Ans3(2,1:25) = out(1:25,2);
%      end
%      if countt==496
%          Ans4(1,1:25) = out(1:25,1);
%          Ans4(2,1:25) = out(1:25,2);
     end
% 

% 
% %close(vidObj);



%% Final Figure 
% h = figure(2);
%             set(h, 'Position', [50 50 1024 640])
%             subplot(10,1,[1 2 3 4 5 6 7 8]);
%             plot(pp(:,1), pp(:,2), 'k', 'LineWidth',1);
%             hold on;
%             plot(Ans1(1,:), Ans1(2,:), '-r+', 'LineWidth',1);
%             plot(Ans2(1,:), Ans2(2,:), '--b', 'LineWidth',1);
%             plot(Ans3(1,:), Ans3(2,:), '--g', 'LineWidth',2);
%             plot(Ans4(1,:), Ans4(2,:), '-m*', 'LineWidth',1);
%             %plot(out(:,1), out(:,2), 'LineWidth',1.5);
% 
%             axis([-3*R 3*R -R/10 2*R]);
%             for q = 0:s:60
%             plot([(q)/div (s+q-1)/div], [-0.03 -.03], '-c', 'LineWidth',4);
%             plot([(q+s-1)/div (q+s)/div], [-0.03 -.03], '-r', 'LineWidth',4);
%             plot([-(s+q-1)/div -(q)/div], [-0.03 -.03], '-c', 'LineWidth',4);
%             plot([-(q+s)/div -(s+q-1)/div], [-0.03 -.03], '-r', 'LineWidth',4);
%             end
%             
%            axis off;
%             daspect([1 1 1]);
%             legend('t* = 0','t* = 0.29','t* = 0.48','t* = 1.45','t* = 1.69');
%             set(gcf,'Color',[1,1,1]);
%             hold off;
%             subplot(10,1,[9 10]);
%             hold on;
%             plot(tt(1:countt-1), lca_arr(1:countt-1),'-r', 'LineWidth',2);
%             plot(tt(1:countt-1), rca_arr(1:countt-1),'-g', 'LineWidth',2);
%             plot(tt(1:countt-1), ca_fitted(1:countt-1),'--k', 'LineWidth',2);
%             axis([0 6 theta_evenl-5 theta_oddl+5]);
%             legend('Actual Left','Actual Right', 'Apparent');
%             set(gcf,'Color',[1,1,1]);
%             xlabel('Time', 'FontSize', 12,'Color', [0 0 0]);
%             ylabel('Contact Angle ( ^o )', 'FontSize', 12,'Color', [0 0 0]);
%             grid on;
%             hold off;

toc;