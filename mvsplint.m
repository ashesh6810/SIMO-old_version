% MVSPLINT  Multivalued Spline Interpolation Function
%   out = SPLINT(P,N) interpolates a cubic spline between points
%   represented by rows of p, N represents the resolution of interpolation
%   and out (NX4) is an array whose first two columns represent
%   interpolated x and y, third column represent slopes at each point on
%   interpolated spline and fourth column represent the respective curvaturs


function [out,normal_vec] = mvsplint(p,N)

length = 0;
params = [0];
for i=2:size(p,1)
    length = length + norm(p(i,1:2)-p(i-1,1:2));
    params = [params;length];
 end
params=params/params(size(params,1));
x = p(:,1);
y = p(:,2);

uu=0:(1/(N-1)):1;
splinex=spline(params,x);
spliney=spline(params,y);
k1=splinex;
k2=params;
SS1=[];

%% calculating xdot from spline X %%

%% calculating ydot from spline Y %%
count1=1;
count2=1;
count3=1;
count4=1;
zz=1;
   for ii=0:(1/(N-1)):1
       xx=0;
       while(xx~=1)
        
       if( ii>=spliney.breaks(count1)&&ii<=spliney.breaks(count1+1))
       
       ydot(zz)=polyval([spliney.coefs(count1,1)*3,spliney.coefs(count1,2)*2,spliney.coefs(count1,3)],ii,[],[spliney.breaks(count1),1]);
       xx=1;
       else
           count1=count1+1;
       end
        if( ii>=spliney.breaks(count2)&&ii<=spliney.breaks(count2+1))
     ydoubledot(zz)=polyval([spliney.coefs(count2,1)*6,spliney.coefs(count2,2)*2],ii,[],[spliney.breaks(count2),1]);
       xx=1;
       else
           count2=count2+1;
        end
        if( ii>=splinex.breaks(count3)&&ii<=splinex.breaks(count3+1))
     xdoubledot(zz)=polyval([splinex.coefs(count3,1)*6,splinex.coefs(count3,2)*2],ii,[],[splinex.breaks(count3),1]);
       xx=1;
       else
           count3=count3+1;
        end
        if( ii>=splinex.breaks(count4)&&ii<=splinex.breaks(count4+1))
       xdot(zz)=polyval([splinex.coefs(count4,1)*3,splinex.coefs(count4,2)*2,splinex.coefs(count4,3)],ii,[],[splinex.breaks(count4),1]);
       xx=1;
       else
           count4=count4+1;
       end
       
       end
     zz=zz+1;
   end
 
 %% slopes %%
 SS=ydot./xdot;

 %% calculating curvature %%
  RR=-(ydoubledot.*xdot-ydot.*xdoubledot)./(xdot.^2+ydot.^2).^1.5;
  %% calculating 2nd derivative %%
  dd=-(ydoubledot.*xdot-ydot.*xdoubledot);
  
 %% calculating values from spline %%
XX = ppval(splinex,uu);
YY = ppval(spliney,uu);
slopes = (gradient(YY,XX));
for i=1:N
    if slopes(i)>=0;
        angles(i)=atand(slopes(i));
    else
        angles(i)=180 +atand(slopes(i));
    end

end
        
dslopes = (gradient(slopes,XX));
kappa = dslopes./(1 + slopes.^2).^(3/2); 
out = zeros(N,4);
out(:,1)=XX;
out(:,2)=YY;
out(:,3)=atand(SS);
out(:,4)=RR;
out(:,5)=kappa;
out(:,6)=dd;
out(:,7)=atand(SS);
out(:,8)=-SS./(XX.*((1+(SS).^2)).^.5);
norm1=cross([zeros(N,2) ones(N,1)], [cosd(out(:,3)), sind(out(:,3)) zeros(N,1)]);

%% setting up the counter  %%
counter=[];
tt=1;
for i=1:N-1
    if out(i,3)>out(i+1,3)
        if abs(out(i,3)-out(i+1,3))>=160
            counter(tt)=i;
            tt=tt+1;
        end
    else
        if abs(out(i+1,3)-out(i,3))>=160
            counter(tt)=i;
            tt=tt+1;
        end
    end
    if (out(i,6)*out(i+1,6))<0
        counter(tt)=i;
        tt=tt+1;
    end
end
k1=counter;
%% updating the counter with the end points %%


  counter_new(1)=1;
  fdg=size(k1,2)+2;
  size(counter_new)
  counter_new(fdg)=N;
  nn=2;
  i=1;
  while(i<=size(counter,2))
      counter_new(nn)=counter(i);
      nn=nn+1;
      i=i+1;
  end
  %% plotting the normals %%
    if out(1,6)>0
   if out(1,3)>0;
      
  kk=sign(-out(1,6));
   else
       kk=sign(out(1,6));
   end
    end
    
    if out(1,6)<0
         if out(1,3)>0;
      
  kk=sign(-out(1,6));
   else
       kk=sign(out(1,6));
       end
    
    end
   
  %% initial sign would depend upon the curvture %%
  cc=1;
  jj=1;

  ii=counter_new(jj);
  while(cc<=size(counter_new,2)-1)
      for i=ii:counter_new(jj+1)
      norm1(i,:)=norm1(i,:)*kk;
      end
      
      ii=counter_new(jj+1)+1;
      jj=jj+1;
      kk=kk*(-1);
      cc=cc+1;
  end
  normal_vec=norm1;
