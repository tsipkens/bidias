
% RUN_INVERSIONS_C  Shorter implementation with customizable regularization parameters.
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

%% Initial guess for iterative schemes
x_init = interp2(grid_b.edges{2}',grid_b.edges{1}',...
reshape(full(b)./(A*ones(size(x0))),grid_b.ne),...
grid_x.elements(:,2),grid_x.elements(:,1));
x_init(isnan(x_init)) = 0;
x_init(isinf(x_init)) = 0;
x_init = sparse(max(0,x_init));
err.init = norm(x0-x_init);


%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
lambda_tk0 = 0.419941123497942;
[x_tk0,D_tk0,L_tk0,Gpo_tk0] = invert.tikhonov(...
    Lb*A,Lb*b,lambda_tk0,0,n_x(1),[],'non-neg');
disp('Inversion complete.');
disp(' ');

err.tk0 = norm(x0-x_tk0);


%% Tikhonov (1st) implementation
disp('Performing Tikhonov (1st) regularization...');
lambda_tk1 = 0.935436889902617;
[x_tk1,D_tk1,L_tk1,Gpo_tk1] = invert.tikhonov(...
    Lb*A,Lb*b,lambda_tk1,1,n_x(1),[],'non-neg');
disp('Inversion complete.');
disp(' ');

err.tk1 = norm(x0-x_tk1);


%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
lambda_tk2 = 1.069019204603001;
[x_tk2,D_tk2,L_tk2] = invert.tikhonov(...
    Lb*A,Lb*b,lambda_tk2,2,n_x(1),[],'non-neg');
disp('Inversion complete.');
disp(' ');

err.tk2 = norm(x0-x_tk2);


%% Twomey
%-- Perform Twomey algorithm ----------------------------%
disp('Performing Twomey...');
x_two = invert.twomey(A,b,x_init,500,[],[],1);

disp('Completed Twomey.');
disp(' ');

err.two = norm(x0-x_two);

