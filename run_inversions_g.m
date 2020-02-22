
% RUN_INVERSIONS_G  Optimize exponential rotated and Tikhonov regularization.
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

%-- Tikhonov (1st order) -----%
disp('Performing Tikhonov (1st) regularization...');
tic;
[x_tk1,lambda_tk1,out_tk1] = optimize.tikhonov_op(...
    Lb*A,Lb*b,[1e-2,1e1],1,n_x(1),x0);
t.tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.tk1 = norm(x0-x_tk1);



%-- Expectation maximization ----%
% Get initial guess
b_init = b;
b_init(b_init<(1e-5*max(b_init))) = 0;
x_init = interp2(grid_b.edges{2}',grid_b.edges{1}',...
    reshape(full(b_init)./(A*ones(size(x0))),grid_b.ne),...
    grid_x.elements(:,2),grid_x.elements(:,1));
x_init(isnan(x_init)) = 0;
x_init(isinf(x_init)) = 0;
x_init = sparse(max(0,x_init));
chi.init = norm(x0-x_init);
x_init_m = grid_x.marginalize(x_init);

disp('Performing expectation-maximization...');
x_em = invert.em(Lb*A,Lb*b,x_init,4,x0);
disp('Inversion complete.');
disp(' ');

chi.em = norm(x0-x_em);
