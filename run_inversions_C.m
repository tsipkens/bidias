
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
chi.init = norm(x0-x_init);


%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
lambda_Tk0 = 0.4199;
[x_Tk0,D_Tk0,L_Tk0,Gpo_Tk0] = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk0,0);
disp('Inversion complete.');
disp(' ');

chi.Tk0 = norm(x0-x_Tk0);


%% Tikhonov (1st) implementation
disp('Performing Tikhonov (1st) regularization...');
lambda_Tk1 = 1.0690;
[x_Tk1,D_Tk1,L_Tk1,Gpo_Tk1] = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk1,1);
disp('Inversion complete.');
disp(' ');

chi.Tk1 = norm(x0-x_Tk1);


%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
lambda_Tk2 = 1.3961;
[x_Tk2,D_Tk2,L_Tk2] = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk2,2);
disp('Inversion complete.');
disp(' ');

chi.Tk2 = norm(x0-x_Tk2);


%% Twomey
%-- Perform Twomey algorithm ----------------------------%
disp('Performing Twomey...');
x_Two = invert.twomey(A,b,x_init,500,[],[],1);

disp('Completed Twomey.');
disp(' ');

chi.Two = norm(x0-x_Two);

