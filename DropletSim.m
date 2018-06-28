clear all;
close all;

% plt = 1 if you want to see the plots as the code is ruuning
plt = 0;
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
    % Deciding the perturbation factor for moving the end points
    perturbation_factor = 1000*N/peri;
    

    
    % rightmost point is on even strip//1 added for symmetry//always negative
    if mod(floor(div*out(1,1)+1), s) == 0 
        delta_R = rca - theta_evenr;
    else
        delta_R = rca - theta_oddr;
    end
    perturbation_R =( R^2*(pi*rca/90 -sind(2*rca))/(perturbation_factor*sind(rca)));
    
    % leftmost point is on even strip
    if mod(floor(div*out(N,1)), s)==0
        delta_L = lca - theta_evenl;
    else
        delta_L = lca - theta_oddl;
    end
    perturbation_L = (R^2*(pi*lca/90 -sind(2*lca))/(perturbation_factor*sind(lca)));

    % decide which side to move
    if delta_R>=0 && delta_L>=0
        outnew = stretch_droplet_sync(out, N, (perturbation_L+perturbation_R)/2, area_orig, 'both',normal_vec);
%     end     % if CA<CA_surface, pin
%     if delta_R>=0 && delta_L<0
%                 outnew = stretch_droplet(out, N, perturbation_R, area_orig, 'right',normal_vec);
%     end
%     if delta_L>=0 && delta_R<0
%                 outnew = stretch_droplet(out, N, perturbation_L, area_orig, 'left',normal_vec);
%     end
%     if delta_L<0 && delta_R<0
    else
        outnew = pin_droplet_sync(out, N, area_orig + Adot*dt,normal_vec);
    end
    [out normal_vec] = mvsplint(outnew,N);

    % calculate the current areas
    area_now  = polyarea(out(:,1),out(:,2));
    area_now_arr(countt)=polyarea(out(:,1),out(:,2));
    % calculate time spent
    t(count) = (area_now-area_first)/Adot;
    
    %counter=counter+1;

    if count>1
    if t(count)>t(count-1)
        % fitting a circle to the shape
        [xc,yc, Re] = circfit(out(:,1),out(:,2));
        xfit = Re*cosd(0:.1:180)+xc; yfit = Re*sind(0:.1:180)+yc;
        ydash = sqrt(Re^2-yc^2)/yc;
        %ca_fitted(countt) = abs(atand(ydash));
        if out(N,3)<=0
            lca_arr(countt)= 180+(out(N,3));
        else
            lca_arr(countt)=out(N,3);
        end
        if atand((ydash))>=0
            ca_fitted(countt)= 180-(atand((ydash)));
        else
            ca_fitted(countt)= abs(atand((ydash)));
        end
        if out(1,3)>=0;
            rca_arr(countt) = 180-out(1,3);
         else
            rca_arr(countt) =abs(out(1,3));                        %abs to eliminate negative - will not work for theta>90
        end
         lca = lca_arr(countt);
         rca = rca_arr(countt);
        tt(countt) = t(count);
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
%    F = getframe(h);
%     writeVideo(vidObj,F);
    end
    end
     if countt==73
         Ans1(1,1:25) = out(1:25,1);
         Ans1(2,1:25) = out(1:25,2);
     end
     if countt==229
         Ans2(1,1:25) = out(1:25,1);
         Ans2(2,1:25) = out(1:25,2);
     end
     if countt==251
         Ans3(1,1:25) = out(1:25,1);
         Ans3(2,1:25) = out(1:25,2);
     end
     if countt==496
         Ans4(1,1:25) = out(1:25,1);
         Ans4(2,1:25) = out(1:25,2);
     end

end

%close(vidObj);



% Final Figure 
h = figure(2);
            set(h, 'Position', [50 50 1024 640])
            subplot(10,1,[1 2 3 4 5 6 7 8]);
            plot(pp(:,1), pp(:,2), 'k', 'LineWidth',1);
            hold on;
            plot(Ans1(1,:), Ans1(2,:), '-r+', 'LineWidth',1);
            plot(Ans2(1,:), Ans2(2,:), '--b', 'LineWidth',1);
            plot(Ans3(1,:), Ans3(2,:), '--g', 'LineWidth',2);
            plot(Ans4(1,:), Ans4(2,:), '-m*', 'LineWidth',1);
%             plot(out(:,1), out(:,2), 'LineWidth',1.5);

            axis([-3*R 3*R -R/10 2*R]);
            for q = 0:s:60
            plot([(q)/div (s+q-1)/div], [-0.03 -.03], '-c', 'LineWidth',4);
            plot([(q+s-1)/div (q+s)/div], [-0.03 -.03], '-k', 'LineWidth',4);
            plot([-(s+q-1)/div -(q)/div], [-0.03 -.03], '-c', 'LineWidth',4);
            plot([-(q+s)/div -(s+q-1)/div], [-0.03 -.03], '-k', 'LineWidth',4);
            end
            
           axis off;
            daspect([1 1 1]);
            legend('t* = 0','t* = 0.29','t* = 0.48','t* = 1.45','t* = 1.69');
            set(gcf,'Color',[1,1,1]);
            hold off;
            subplot(10,1,[9 10]);
            hold on;
            plot(area_now_arr(1:countt-1)/pitch^2, lca_arr(1:countt-1),'-r', 'LineWidth',2);
            plot(area_now_arr(1:countt-1)/pitch^2, rca_arr(1:countt-1),'-g', 'LineWidth',2);
            plot(area_now_arr(1:countt-1)/pitch^2, ca_fitted(1:countt-1),'--k', 'LineWidth',2);
            axis([0 300 theta_evenl-5 theta_oddl+5]);
            legend('Actual Left','Actual Right', 'Apparent');
            set(gcf,'Color',[1,1,1]);
            xlabel('Time', 'FontSize', 12,'Color', [0 0 0]);
            ylabel('Contact Angle ( ^o )', 'FontSize', 12,'Color', [0 0 0]);
            grid on;
            hold off;

toc;
