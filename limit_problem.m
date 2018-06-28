%% plots for comparisons %%

clear all;
close all;
NN=30;
N=100;

theta=110;
R=2.5;
%% control points %%
alpha=90-theta:2*theta/NN:90+theta;

pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];
%% passing it through mvsplint %%    
  [out,coefsX,coefsY,gradX,gradY,U,kk,splinex,spliney]=splinetest(pp,N);
  figure(1)
  plot(out(:,1),out(:,2));
  daspect([1 1 1])
  ylabel('shape')
  figure(2)
  plot(U,out(:,4));
   ylabel('curvature from mvsplint')
  figure(10)
  plot(U,out(:,3));
    ylabel('slope from mvsplint')
 
  
 %% evaluating xdoubledot%%
  count=1;
 zz=1;
 for ii=0:(1/(N-1)):1
     xx=0;
      while(xx~=1)
       if( ii>=splinex.breaks(count)&&ii<=splinex.breaks(count+1))
     xdoubledot(zz)=polyval([splinex.coefs(count,1)*6,splinex.coefs(count,2)*2],ii,[],[splinex.breaks(count),1]);
       xx=1;
       else
           count=count+1;
       end
       
       end
  zz=zz+1;  
 end
 figure(3);
 plot(U,xdoubledot);
 ylabel('xdot')
%% evaluating ydoubledot %%
 count=1;
 zz=1;
 for ii=0:(1/(N-1)):1
     xx=0;
      while(xx~=1)
       if( ii>=spliney.breaks(count)&&ii<=spliney.breaks(count+1))
     ydoubledot(zz)=polyval([spliney.coefs(count,1)*6,spliney.coefs(count,2)*2],ii,[],[spliney.breaks(count),1]);
       xx=1;
       else
           count=count+1;
       end
       
       end
  zz=zz+1;  
 end
 figure(4);
 plot(U,ydoubledot);
 ylabel('ydot')
 %% evaluating xdot %%
count=1;
zz=1;
   for ii=0:(1/(N-1)):1
       xx=0;
       while(xx~=1)
        
       if( ii>=splinex.breaks(count)&&ii<=splinex.breaks(count+1))
       
       xdot(zz)=polyval([splinex.coefs(count,1)*3,splinex.coefs(count,2)*2,splinex.coefs(count,3)],ii,[],[splinex.breaks(count),1]);
       xx=1;
       else
           count=count+1;
       end
       end
     zz=zz+1;
   end
 figure(5);
 plot(U,xdot);
  ylabel('xdoubeldot')
%% evaluating ydot %%
count=1;
zz=1;
   for ii=0:(1/(N-1)):1
       xx=0;
       while(xx~=1)
        
       if( ii>=spliney.breaks(count)&&ii<=spliney.breaks(count+1))
       
       ydot(zz)=polyval([spliney.coefs(count,1)*3,spliney.coefs(count,2)*2,spliney.coefs(count,3)],ii,[],[spliney.breaks(count),1]);
       xx=1;
       else
           count=count+1;
       end
       end
     zz=zz+1;
   end
 figure(6);
 plot(U,ydot);
 ylabel('ydoubeldot')
 %% evaluating slope from polynomial %%
 slope=ydot./xdot;
 figure(7)
 plot(U,slope)
 ylabel('slope from polynomial')
 %% evaluating double derivative from polynomial %%
 dslope=(ydoubledot.*xdot-ydot.*xdoubledot);
 figure(8)
 plot(U,dslope);
  ylabel('2nd derivative from polynomial')
  %% evaluating curvature %%
  RR=-(ydoubledot.*xdot-ydot.*xdoubledot)./(xdot.^2+ydot.^2).^1.5;
  figure(9)
  plot(U,RR);
  ylabel('curvature from polynomial')
  %%
%     figure(3);
%     plot(U,out(:,4));
%       ylabel('curvature from mvsplint')
%     figure(11);
%     plot(out(:,1),out(:,2));
%     daspect([1 1 1])
%     ylabel('shape')
%     xx=uu_upper:.01:uu_lower
%     xdot=polyval(c1,xx,[],[uu_upper,uu_lower]);
%     figure(1)
%     plot(xx,xdot);
%     ylabel('xdot')
%     cc1(1)=coefsY(N,1)*3;
%     cc1(2)=coefsY(N,2)*2;
%     cc1(3)=coefsY(N,3)*1;
%     yy=uu_upper:.01:uu_lower
%      ydot=polyval(cc1,yy,[],[uu_upper,uu_lower]);
%      figure(2)
%      plot(yy,ydot)
%        ylabel('ydot')
%      Z=ydot./xdot;
%      figure(4)
%      plot(yy,Z)
%        ylabel('slope from polynomial')
%      c2(1)=coefsX(N,1)*6;
%      c2(2)=coefsX(N,2)*2;
%      xdoubledot=polyval(c2,xx,[],[uu_upper,uu_lower]);
%      figure(5)
%      plot(xx,xdoubledot)
%        ylabel('xdoubledot')
%        cc2(1)=coefsY(N,1)*6;
%        cc2(2)=coefsY(N,2)*2;
%        ydoubledot=polyval(cc2,yy,[],[uu_upper,uu_lower]);
%        figure(6)
%        plot(yy,ydoubledot)
%          ylabel('ydoubledot')
%        ZD=ydoubledot./xdot.^2 -(ydot.*xdoubledot)/xdot.^3;
%        figure(7)
%        plot(xx,ZD)
%          ylabel('2nd derivative with polynomial')
%        RR=(ydoubledot.*xdot-ydot.*xdoubledot)./(xdot.^2+ydot.^2).^1.5;
%        figure(8)
%        plot(xx,RR)
%          ylabel('curvature from polynomial')
%          figure(12)
%          plot(U,out(:,3))
%        