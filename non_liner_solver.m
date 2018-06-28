function angle = non_liner_solver(out,N,area,lca)
ashesh=1;


while(ashesh==1||( d-area)>.001)
    ashesh=2;
    d=((pi/180)*lca*(cscd(lca))^2-cotd(lca))*(out(1,1)-(out(1,1)+out(N,1))/2)^2;
   lca=lca-.0001;
end
angle=lca;
end