% MVSPLINT  Multivalued Spline Interpolation Function
%   out = SPLINT(P,N) interpolates a cubic spline between points
%   represented by rows of p, N represents the resolution of interpolation
%   and out (NX4) is an array whose first two columns represent
%   interpolated x and y, third column represent slopes at each point on
%   interpolated spline and fourth column represent the respective curvaturs


function out = drthakur(p,N)

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

XX = ppval(splinex,uu);
YY = ppval(spliney,uu);
slopes = gradient(YY,XX);
dslopes = (gradient(slopes,XX));
kappa = dslopes./(1 + slopes.^2).^(3/2);


 
out = zeros(N,4);
out(:,1)=XX;
out(:,2)=YY;
out(:,3)=slopes;
out(:,4)=kappa;

