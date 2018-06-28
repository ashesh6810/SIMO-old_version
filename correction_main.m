%% correction main %%
% clear all;
% NN=30;
% theta=60; %% any angle of choice for a circle %%
% R=2.5;
% alpha=90-theta:2*theta/NN:90+theta;
% 
% pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];
% [out,U]=mvsplint(pp,NN);
% norm=cross([zeros(length(U),2) ones(length(U),1)], [cosd(out(:,3)), sind(out(:,3)) zeros(length(U),1)]);
% 
% figure(3)
% plot(out(:,1),out(:,2));
% figure(1)
% quiver(out(:,1),out(:,2),norm(:,1),norm(:,2));
% hold on

function[normal_vec]= correction_main(out,U,NN)
%% calculating the initial normals %%
norm=cross([zeros(length(U),2) ones(length(U),1)], [cosd(out(:,3)), sind(out(:,3)) zeros(length(U),1)]);

%% setting up the counter  %%
counter=[];
tt=1;
for i=1:length(U)-1
    if out(i,3)>out(i+1,3)
        if abs(out(i,3)-out(i+1,3))>=160
            counter(tt)=i;
            tt=tt+1;
        end
    else
        if abs(out(i+1,3)-out(i,3))>=160
            counter(tt)=i;
        end
    end
    if (out(i,6)*out(i+1,6))<0
        counter(tt)=i;
        tt=tt+1;
    end
end
%% updating the counter with the end points %%

 
  counter_new(1)=1;
  counter_new(length(counter)+2)=NN;
  nn=2;
  i=1;
  while(i<=length(counter))
      counter_new(nn)=counter(i);
      nn=nn+1;
      i=i+1;
  end
  %% plotting the normals %%
   jj=1;
  kk=1; %% initial sign would depend upon the curvture %%
  cc=1;
  m=1;
  ii=counter_new(jj);
  while(cc<=length(counter_new)-1)
      for i=ii:counter_new(jj+1)
      norm(i,:)=norm(i,:)*kk;
      end
      
      ii=counter_new(jj+1)+1;
      jj=jj+1;
      kk=kk*(-1);
      cc=cc+1;
  end
  normal_vec=norm;
%         figure(2)
%         quiver(out(:,1),out(:,2),norm(:,1),norm(:,2));
        
