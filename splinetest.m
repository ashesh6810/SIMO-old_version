%% testspline %%
function[ out,k1,k2,k3,k4,k5,k6,k7,k8] = splinetest(p,N)

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

k5=uu;
k7=splinex;
k8=spliney;
XX = ppval(splinex,uu);
YY = ppval(spliney,uu);
slopes = gradient(YY,XX);
k3=gradient(XX,uu);
k4=gradient(YY,uu);

dslopes = (gradient(slopes,XX));
kappa = dslopes./(1 + slopes.^2).^(3/2);
out = zeros(N,5);
out(:,1)=XX;
out(:,2)=YY;
out(:,3)=slopes;
out(:,4)=kappa;
out(:,5)=dslopes;
k1=splinex.coefs;
k2=spliney.coefs;
k6=params;