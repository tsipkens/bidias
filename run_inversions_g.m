
% RUN_INVERSIONS_G  Optimize exponential distance and 1st-order Tikhonov regularization.
% Author: Timothy Sipkens, 2019-05-28
%=========================================================================%

%-- Tikhonov (1st order) -----%
disp('Performing Tikhonov (1st) regularization...');
tic;
[x_tk1,lambda_tk1,out_tk1] = optimize.tikhonov_op(...
    Lb*A,Lb*b,[1e-2,1e1],1,n_x(1),x0);
t.tk1 = toc;
disp('Inversion complete.');
disp(' ');

eps.tk1 = norm(x0-x_tk1);


%%
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
eps.init = norm(x0-x_init);
x_init_m = grid_x.marginalize(x_init);

disp('Performing expectation-maximization (attempt no. 1)...');
x_em = invert.em(Lb*A,Lb*b,x_init,3,x0);
disp('Inversion complete.');
disp(' ');

eps.em = norm(x0-x_em);



%-- Expectation maximization (attempt 2) ----%
% Get initial guess
n2 = [18,21];
grid2 = Grid([grid_t.span],...
    n2,'logarithmic');
B2 = grid2.transform(grid_t); % evaluate matrix modifier to transform kernel
A2 = A_t*B2;
x02 = grid2.project(grid_t,x_t); % project into basis for x

x_init2 = interp2(grid_b.edges{2}',grid_b.edges{1}',...
    reshape(full(b_init)./(A2*ones(size(x02))),grid_b.ne),...
    grid2.elements(:,2),grid2.elements(:,1));
x_init2(isnan(x_init)) = 0;
x_init2(isinf(x_init2)) = 0;
x_init2 = sparse(max(0,x_init2));
eps.init2 = norm(x02-x_init2);

disp('Performing expectation-maximization (attempt no. 2)...');
x_em2 = invert.em(Lb*A2,Lb*b,x_init2,250);
disp('Inversion complete.');
disp(' ');

eps.em2 = norm(x02-x_em2);




%-- Expectation maximization (attempt 3) ----%
% Get initial guess
x_init3 = ones(size(x02));

disp('Performing expectation-maximization (attempt no. 3)...');
x_em3 = invert.em(Lb*A2,Lb*b,x_init3,250);
disp('Inversion complete.');
disp(' ');

eps.em3 = norm(x02-x_em3);



%-- Expectation maximization (attempt 4) ----%
disp('Performing expectation-maximization (attempt no. 1)...');
x_em4 = invert.em(Lb*A,Lb*b,ones(size(x_init)),3,x0);
disp('Inversion complete.');
disp(' ');

eps.em4 = norm(x0-x_em4);



%-- Plot EM results ------------------------%
%{
figure(31);
grid_x.plot2d(x_em);
colormap(cm);
colorbar;
title('EM, attempt 1 (x_{init} = interp(b_{i}/a_{i}x)');

figure(32);
grid2.plot2d(x_em2);
colormap(cm);
colorbar;
title('EM, attempt 2 (x_{init} = interp(b_{i}/a_{i}x)');

figure(33);
grid2.plot2d(x_em3);
colormap(cm);
colorbar;
title('EM, attempt 3 (x_{init} = 1)');
%}





