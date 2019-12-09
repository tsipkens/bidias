
% RUN_INVERSIONS_I  Multi-dimensional optimized exponential rotated and exp. dist. regularization
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

guess = [1.3,1/4,log10(1.8),0.84]; % [lambda, ratio, ld, corr]
disp('Optimizing exponential distance regularization...');
[x_exp_opt,lambda_exp_opt,out_exp_opt] = optimize.exp_dist_opx(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    guess,x0); 
disp('Inversion complete.');
disp(' ');


disp('Parametric study of exponential distance regularization...');
[x_exp_par,lambda_exp_par,out_exp_par] = optimize.exp_dist_opbf(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    guess,x0); 
disp('Inversion complete.');
disp(' ');

chi.exp_dist = norm(x_exp_opt-x0);


%%
%{
[vec_d1,vec_d2] = ndgrid(grid_x.elements(:,2),grid_x.elements(:,2));
[vec_m1,vec_m2] = ndgrid(grid_x.elements(:,1),grid_x.elements(:,1));
Gd_fun = @(y) [(y(3)/y(2))^2,min(y(4)*y(3)^2/y(2),0.98);...
    min(y(4)*y(3)^2/y(2),0.98),y(3)^2]; % version for no correlation
    % y(2) = ratio, y(3) = ld, y(4) = corr

out = out_exp_par;
fprintf('Evaluating F:');
tools.textbar(0);
for kk=1:length(out)
    y = [out(kk).lambda,...
        out(kk).ratio,...
        out(kk).ld,...
        out(kk).corr_vec];
    
    Gd = Gd_fun(y);
    Gd_inv = inv(Gd);
    drm = log10(vec_m1)-log10(vec_m2);
    drd = log10(vec_d1)-log10(vec_d2);
    d = sqrt(drm.^2.*Gd_inv(1,1)+...
        2.*drd.*drm.*Gd_inv(1,2)+...
        drd.^2.*Gd_inv(2,2)); % distance
    
    Gpr = exp(-d);
    
    Gpr_inv = pinv(Gpr);
    [Lpr,~] = chol(Gpr_inv);
    clear Gpr_inv; % to save memory
    Lpr = y(1).*Lpr./max(max(Lpr));
    % Lpr(abs(Lpr)<(0.01.*mean(mean(abs(Lpr))))) = 0;
    Lpr(d>1.5) = 0;
    Lpr = sparse(Lpr);
    
    fit_b(kk) = norm(Lb*(A*out(kk).x-b));
    fit_pr(kk) = norm(Lpr*out(kk).x);
    
    tools.textbar(kk/length(out_exp_par))
end
% B = 1/2.*(-(fit_pr+fit_b) -det_pr+det_po);
F = -1/2.*(fit_pr+fit_b);
%}
