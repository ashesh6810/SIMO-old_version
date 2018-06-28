function outnew = stretch_droplet_Bo_axis(out, N, perturbation, area_orig, side,kk,Bo)
% Decise which end to move and then move it
if strcmp(side,'right')
    out(1,1) = out(1,1) + perturbation/norm(out(2,1:2)-out(1,1:2));    
end
if strcmp(side,'left')
    out(N,1) = out(N,1) - perturbation/norm(out(N,1:2)-out(N-1,1:2));    
end
if strcmp(side,'both')
    out(1,1) = out(1,1) + perturbation/norm(out(2,1:2)-out(1,1:2));   
    out(N,1) = out(N,1) - perturbation/norm(out(N,1:2)-out(N-1,1:2));    
end

% Calculate the current area
%area_curr = polyarea(out(:,1),out(:,2));
outvolume=mvsplint(out,1000);
area_curr=pappus(outvolume);
% Calculate the dA
del_area = area_orig - area_curr;

% Set the PID Parameters
P = 10*sqrt(abs(del_area)/2);
I = P/100;
D = P/10;

%Initialize
err_int = 0;
err_prev = 0;
err_tol = area_orig*1e-5;

%% Run the loop for moving anyother point such that dA due to end point movement is 0
while abs(del_area) > err_tol  
    % calculate the surface normal

    % multiply the surface normal with magnitude and the direction of curvature 

      
for i=1:N
kk(i,:)=abs(out(i,4)).*kk(i,:);
end

    k1=kk;
    comp=(out(2:N-1,4)+out(2:N-1,8))-Bo*out(2:N-1,2);
    % find the inex of the point with minimum curvature
    [minimum, min_ind] = min(comp);
    % end points are constrained see a line before
    min_ind = min_ind + 1;
    % decide the componenents in each direction
    perturb_vec = kk(min_ind, 1:2);
      perturb_vec_sync = [-kk(min_ind,1),kk(min_ind,2)];
    % normalize these
    perturb_vec = perturb_vec/norm(perturb_vec);
    perturb_vec_sync = perturb_vec_sync/norm(perturb_vec_sync);
    % move this point using the PID controller
    out(min_ind,1:2) = out(min_ind,1:2) - P*(del_area)/2*perturb_vec - ...
        I*(err_int)/2*perturb_vec...   
        -D*(del_area - err_prev)/2*perturb_vec;  
    
     out(N-min_ind+1,1:2) = out(N-min_ind+1,1:2) - P*(del_area)/2*perturb_vec_sync - ...
        I*(err_int)/2*perturb_vec_sync...   
        -D*(del_area - err_prev)/2*perturb_vec_sync;  
    
    
    
    
    % store the previous dA
    err_prev = del_area;
    % update the points with a new spline
   [ out_updated,normal_vec] = mvsplint(out,N);    
    out = out_updated;
    outvolume=mvsplint(out,1000);
    
    kk=normal_vec;
%     k2=k1;
    % calculate the current area
   % area_curr=polyarea(out_updated(:,1),out_updated(:,2));
   area_curr=pappus(outvolume);
    % calculate the new dA
    del_area = area_orig - area_curr;
    % compute the error integral term
    err_int = err_int + (del_area);
    % initialize integral term if the error is zero
    if err_int>2*area_curr
        err_int = 0;
    end
    
end

% return the out matrix
outnew = out;

end