function outnew = pin_droplet(out, N, area_orig,normal_vec)

% Calculate the current area
area_curr = polyarea(out(:,1),out(:,2));

% Calculate the dA
del_area = area_orig - area_curr;

% Set the PID Parameters
P = 10*sqrt(abs(del_area));
I = P/100;
D = P/10;

%Initialize
err_int = 0;
err_prev = 0;
err_tol = area_orig*1e-5;

%% Run the loop until the new area accomodates dA, volume addition
while abs(del_area) > err_tol

    for i=1:N
        normal_vec(i,:) = (out(i,4)).*normal_vec(i,:);
    end
        
    % find the inex of the point with maximum curvature
    [maximum, max_ind] = max(out(2:N-1,4));
    % end points are constrained see a line before
    max_ind = max_ind + 1;
    % decide the componenents in each direction
    perturb_vec = normal_vec(max_ind, 1:2);
    % normalize these
    perturb_vec = perturb_vec./norm(perturb_vec);
    % move this point using the PID controller
    out(max_ind,1:2) = out(max_ind,1:2) + P*(del_area)*perturb_vec + ...
        +I*(err_int)*perturb_vec...
        +D*(del_area - err_prev)*perturb_vec;
    % store the previous dA
    err_prev = del_area;
    % update the points with a new spline
    [out_updated,normal_vec] = mvsplint(out,N);    
    out = out_updated;
    % calculate the current area
    area_curr=polyarea(out_updated(:,1),out_updated(:,2));
    % calculate the new dA
    del_area = area_orig - area_curr;
    % compute the error integral term
    err_int = err_int + (del_area);
    % initialize integral term if the error is zero
    if err_int>area_curr
        err_int = 0;
    end
    
end

% return the out matrix
outnew = out;



