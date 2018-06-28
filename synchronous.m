clear all
close all
%% initial coordinates %%
X=[0,1,2,3,4,5,6,7,8,9,10]./1.5;
X1=[1,2,3,4,5,6,7,8,9]/1.5;
%% free to move coordinates %%
Y=[ 0 1.5 4 6 5 3 5 6 4 1.5 0]./1.5;
Y1=[1.5 4 6 5 3 5 6 4 1.5]./1.5;
%% definig the p matrcies for dr thakurs function %%
pp(:,1)=X;
pp(:,2)=Y;
%% taking N points for runiing the algo and H points to plot the shape %%
N=60;
H=100;
epsilon=.01
out=drthakur(pp,N);
UU=drthakur(pp,H);
%% computing the initial area of the shape %%
Z=polyarea(UU(:,1),UU(:,2));
%% calculating initial perimeter %%
p4=1;
for i=1:N-1
    P1(p4)=[sqrt(1+(out(i,3))^2)]*(abs(out(i+1,1))-abs(out(i,1)));
    p4=p4+1;
end
jj=1;
perimeter=0;

for i=1:N-1
    perimeter=perimeter+P1(i);
end

%% applying the algoirthm for N-2 points %%
for i=2:N-1
    out1(i-1,1)=out(i,1);
    out1(i-1,2)=out(i,2);
    out1(i-1,3)=out(i,3);
    out1(i-1,4)=out(i,4);
end
%% defining a separate R matrix for the curvature %%
p3=1;
for i=1:N-2
    R(p3)=out1(i,4);
    p3=p3+1;
end
%% finding out the indices for the maximum and minimum value of curvature %%
[l,mx]=max(out1(:,4));
[n,mn]=min(out1(:,4));
Rmax=[];
Rmin=[];
Rmax_ind=[];
Rmin_ind=[];
cnt=1;
cntt=1;
for i=1:N-2
    
    if abs(l-out1(i,4))<(l-n)/1000
       Rmax(i)=out(i,4);
       Rmax_ind(cnt)=i;
       cnt=cnt+1;
    end
    
    if abs(n-out1(i,4))<(l-n)/1000
        Rmin(i)=out(i,4);
        Rmin_ind(cntt)=i;
        cntt=cntt+1;
    end
end

%% starting the algorithm %%;
%% defining initial parameters
alpha=.002;
count=1;
ashesh=1;
%%

     while(ashesh==1||(perimeter1(count-100)>perimeter1(count-1)))
%          for M=1:2000
        if count<=100
         ashesh=1;
        else
    ashesh=2;
        end
    si_y_min=sign(R(mn));                                % compressing the 4 conditions for X %
    si_x_min=-sign(out1(mn,3))/sign(R(mn));              % compressing the 4 conditions for Y %
    
    si_y_max=sign(R(mx));                                % compressing the 4 conditions for X %
    si_x_max=-sign(out1(mx,3))/sign(R(mx));              % compressing the 4 conditions for y %
    
    for i=1:length(Rmin_ind)
    out1(Rmin_ind(i),1)=out1(Rmin_ind(i),1)+(si_x_min)*abs(out1(Rmin_ind(i),3))/sqrt(1+out1(Rmin_ind(i),3)^2)*alpha;
    out1(Rmin_ind(i),2)=out1(Rmin_ind(i),2)+((si_y_min))*1/sqrt(1+out1(Rmin_ind(i),3)^2)*alpha;
    end
    for i=2:N-1
        out(i,2)=out1(i-1,2);
        out(i,2)=out1(i-1,2);
    end
    
    AA=polyarea(out(:,1),out(:,2));
    
    k=.01;
    ash=1;
    bb=1;
   
    while(ash==1||(abs(Z-b)> Z/1000))
        ash=2;
        for i=1:length(Rmax_ind)
        out1(Rmax_ind(i),1)=out1(Rmax_ind(i),1)-(out1(Rmax_ind(i),3))/sqrt(1+out1(Rmax_ind(i),3)^2)*k;
        out1(Rmax_ind(i),2)=out1(Rmax_ind(i),2)+1/sqrt((out1(Rmax_ind(i),3))^2+1)*k;
        end
        for i=2:N-1
            out(i,1)=out1(i-1,1);
        end
        for i=2:N-1
            out(i,2)=out1(i-1,2);
        end
        
        
        b= polyarea(out(:,1),out(:,2));
%         if M==557 
%         BB(bb)=polyarea(out(:,1),out(:,2));
%     end
        if abs(Z-b)> Z/1000
            for i=1:length(Rmax_ind)
            out1(Rmax_ind(i),1)=out1(Rmax_ind(i),1)+(out1(Rmax_ind(i),3))/sqrt(1+(out1(Rmax_ind(i),3))^2)*k;
            out1(Rmax_ind(i),2)=out1(Rmax_ind(i),2)-1/sqrt((out1(Rmax_ind(i),3))^2+1)*k;
            end
                k=k+(Z-b)/1000
           
            
        end
        
        
        
        bb=bb+1;
    end
    ash=3;
    
    %% algorithm ends %%
    
    G=[];
    for i=1:N
        G(i,1)=out(i,1);
    end
    for i=1:N
        G(i,2)=out(i,2);
    end
    %% plotting the points %%
    out=drthakur(G,N);
    UU=drthakur(G,H);
    f=spapi(3,UU(:,1),UU(:,2));
    figure(3)
    subplot(2,1,2)
%     fnplt(f,2);
plot(UU(:,1),UU(:,2));
    axis([0 10,0 10]);
    axis equal
    hold off
    %% plotting the curvatures %%
    for i=1:N-1
        D(i)=[sqrt(1+(out(i,3))^2)]*(abs(out(i+1,1))-abs(out(i,1)));
    end
    U=[];
    U(1)=0;

    for i=2:N
        U(i)=U(i-1)+D(i-1);
%         J(i)=abs((.4403+7674)/2-U(i));
    end
    for i=1:N
        RR(i)=out(i,4);
    end
    h=spapi(2,U,RR);
    
        ff=1;
        for i=N:-1:1
            vv(ff)=RR(i);
            ff=ff+1;
        end
%         jj=spapi(2,uu,vv);
    figure(3)
%     if M==1
%         fnplt(h,2)
%         axis([0 15 -7.5 7.5])
%         hold on
%     else
subplot(2,1,1)
        plot(U,RR,':b*',U,vv,':r*');
        hold on;
     
        axis([0 15 -7.5 7.5])
        hold off
        
        

     
    %% re-iterating the floating points and the optimum indices %%
    for i=2:N-1
        out1(i-1,1)=out(i,1);
        out1(i-1,2)=out(i,2);
        out1(i-1,3)=out(i,3);
        out1(i-1,4)=out(i,4);
    end
    
    p3=1;
    for i=2:N-1
        R(p3)=out(i,4);
        p3=p3+1;
    end
    [l,mx]=max(out1(:,4));
    [n,mn]=min(out1(:,4));
    Rmax=[];
Rmin=[];
Rmax_ind=[];
Rmin_ind=[];
cnt=1;
cntt=1;
for i=1:N-2
    
    if abs(l-out1(i,4))<(l-n)/1000
       Rmax(i)=out(i,4);
       Rmax_ind(cnt)=i;
       cnt=cnt+1;
    end
    
    if abs(n-out1(i,4))<(l-n)/1000
        Rmin(i)=out(i,4);
        Rmin_ind(cntt)=i;
        cntt=cntt+1;
    end
end
    
    p4=1;
    for i=1:N-1
        P1(p4)=[sqrt(1+(out(i,3))^2)]*(abs(out(i+1,1))-abs(out(i,1)));
        p4=p4+1;
    end
    perimeter1(count)=0;
    
    for i=1:N-1
        perimeter1(count)=perimeter1(count)+P1(i);
    end
    count=count+1;
    
    
    %% to be deleted %%
     if count==8
    figure(1)
    subplot(2,1,1)
    plot(U,RR,':b*');
   axis([ 0 16 -5 5])
   set(gca,'DataAspectRatio',[1 1 1])
    hold on
%     figure(2)
    subplot(2,1,2)
    fnplt(f,2,'b')
    hold on
       axis([0 10 0 10])
    set(gca,'DataAspectRatio',[1 1 1])
    end
    if count==500
    figure(1)
    subplot(2,1,1)
    plot(U,RR,':b');
     axis([ 0 16 -5 5])
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
%     figure(2)
    subplot(2,1,2)
    fnplt(f,2,'k')
    
    
     axis([0 10 0 10])
    hold on
    end
    
    if count==2000
    figure(1)
    subplot(2,1,1)
     plot(U,RR,'b-');
   axis([ 0 16 -5 5])
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
%     figure(2)
    subplot(2,1,2)
    fnplt(f,2,'c')
    hold on
       axis([0 10 0 10])
    set(gca,'DataAspectRatio',[1 1 1])
    end
    
%     if M==290
%     figure(1)
% %     subplot(2,1,1)
%     fnplt(h,2,'r')
%     hold on
%     figure(2)
% %     subplot(2,1,2)
%     fnplt(f,2,'r')
%     hold on
%      axis([0 10 0 10])
%     axis equal
%     end
    
    if count==8900
    figure(1)
    subplot(2,1,1)
      plot(U,RR,'b+')
     axis([ 0 16 -5 5])
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
%     figure(2)
    subplot(2,1,2)
    fnplt(f,2,'y')
    hold on
    axis([0 10 0 10])
    set(gca,'DataAspectRatio',[1 1 1])
    end
     end   %%remember%%
% end
figure(3)
subplot(2,1,2)
fnplt(f,2)
axis([0 10 0 10])
axis equal
hold on
% theta=atan(out(1,3));
%    z1=1;
%    n=100;
% for i=0:10/n:10                    % defining the circle X points %
%     X1(z1)=i;
%     z1=z1+1;
% end
% for i=1:n+1
%     Y1(i)=sqrt(25*(csc(theta))^2-(X1(i)-5)^2)-5*cot(theta);   % defining the y Points%
% end
% g=spapi(3,X1,Y1);
% figure(1)
% fnplt(g,2,'k--');
% axis([0 10 0 10]);
% hold on
%
% axis equal
% pt1=[out(1,1),out(1,2)];
% pt2=[out(30,1),out(30,2)];
% pt3=[out(60,1),out(60,2)];
% [centre,radius]=calcCircle(pt1,pt2,pt3);
% 
% R = radius;
% Center = [centre(1),centre(2)];
% CIRCLE(Center,R,1000,'b-');
% axis([0 10 0 10])
% hold on
% figure(1)
% plot(Center(1),Center(2),'g.')
% axis([0 10 0 10])



XY=[];
XY(:,1)=out(:,1);
XY(:,2)=out(:,2);
PAR=CircleFitByTaubin(XY);
centre(1)=PAR(1);
centre(2)=PAR(2);
RADIUS=PAR(3);
CIRCLE(centre,RADIUS,1000,'b-');
axis([0 10 0 10])
hold on
figure(3)
subplot(2,1,2)
plot(centre(1),centre(2),'r');
axis([0 10 0 10])