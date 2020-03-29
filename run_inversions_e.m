
% RUN_INVERSIONS_E  Run full set of optimal inversions for CPMA-SP2 inversion.
% Author:           Arash Naseri, Timothy Sipkens, 2020-02-06
%=========================================================================%

%% Initial guess for iterative schemes
b_init = b;
b_init(b_init<(1e-5*max(b_init))) = 0;
x_init = interp2(grid_b.edges{2}',grid_b.edges{1}',...
    reshape(full(b_init)./(A*ones(size(x0))),grid_b.ne),...
    grid_x.elements(:,2),grid_x.elements(:,1));
x_init(isnan(x_init)) = 0;
x_init(isinf(x_init)) = 0;
x_init = sparse(max(0,x_init));
eps.init = 100*norm(x0-x_init)/norm(x0);
x_init_m = grid_x.marginalize(x_init);


%% Least squares
disp('Performing LS inversion...');
x_length = length(A(1,:));
x_lsq = invert.lsq(Lb*A,Lb*b);
disp('Inversion complete.');
disp(' ');

eps.lsq = 100*norm(x0-x_lsq)/norm(x0);


%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
tic;
[x_tk0,lambda_tk0,out_tk0] = optimize.tikhonov_op(...
    Lb*A,Lb*b,[1e-2,1e2],0,grid_x,...
    x0,[],'interior-point');
t.tk0 = toc;
disp('Inversion complete.');
disp(' ');

eps.tk0 = 100*norm(x0-x_tk0)/norm(x0);


%% Tikhonov (1st) implementation
disp('Performing Tikhonov (1st) regularization...');
tic;
[x_tk1,lambda_tk1,out_tk1] = optimize.tikhonov_op(...
    Lb*A,Lb*b,[1e-2,1e2],1,grid_x,...
    x0,[],'interior-point');
t.tk1 = toc;
disp('Inversion complete.');
disp(' ');

eps.tk1 = 100*norm(x0-x_tk1)/norm(x0);


%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
tic;
[x_tk2,lambda_tk2,out_tk2] = optimize.tikhonov_op(...
    Lb*A,Lb*b,[1e-2,1e2],2,grid_x,...
    x0,[],'interior-point');
t.tk2 = toc;
disp('Inversion complete.');
disp(' ');

eps.tk2 = 100*norm(x0-x_tk2)/norm(x0);


%% Tikhonov (0th) implementation L-curve
disp('L-curve optimization of Tikhonov (0th) regularization...');
tic;
Lpr0 = invert.tikhonov_lpr(0,grid_x);
[x_tk0_lc,lambda_tk0_lc] = optimize.tikhonov_lcurve(...
    Lb*A,Lb*b,[10^-3,100],Lpr0);
t.tk0_lc = toc;
disp('Inversion complete.');
disp(' ');
eps.tk0_lcurve = 100*norm(x0-x_tk0_lc)/norm(x0);


%% Tikhonov (1st) implementation L-curve
disp('L-curve optimization of Tikhonov (1st) regularization...');
tic;
Lpr0 = invert.tikhonov_lpr(1,grid_x);
[x_tk1_lc,lambda_tk1_lc,res_norm,x_norm] = ...
    optimize.tikhonov_lcurve(...
    Lb*A,Lb*b,[10^-3,100],Lpr0);
t.tk1_lc = toc;
disp('Inversion complete.');
disp(' ');
eps.tk1_lcurve = 100*norm(x0-x_tk1_lc)/norm(x0);


%% Tikhonov (2nd), L-curve
disp('L-curve optimization of Tikhonov (2nd) regularization...');
tic;
Lpr0 = invert.tikhonov_lpr(2,grid_x);
[x_tk2_lc,lambda_tk2_lc] = optimize.tikhonov_lcurve(...
    Lb*A,Lb*b,[10^-3,100],Lpr0);
t.tk2_lc = toc;
disp('Inversion complete.');
disp(' ');
eps.tk2_lc = 100*norm(x0-x_tk2_lc)/norm(x0);


%% MART, Maximum entropy regularized solution

disp('Performing MART...');
tic;
[x_mart,iter_mart,out_mart] = ...
    optimize.mart_op(A,b,x_init,1:300,x0);
t.MART = toc;
disp('Inversion complete.');
disp(' ');

eps.mart = 100*norm(x0-x_mart)/norm(x0);


%% MART with smoothing function
disp('Performing MART with smoothing function...');
tic;
[x_mart_mh,Sf_mart_mh,out_mart_mh] = optimize.martmark_op(...
    A,b,Lb,grid_x,x_init,35,[1e2,1e-5],x0,'Buckley');
t.mart_mh = toc;

eps.mart_mh = 100*norm(x0-x_mart_mh)/norm(x0);

disp('Completed MART with smoothing function.');
disp(' ');


%% Twomey
disp('Performing Twomey...');
tic;
[x_two,iter_two,out_two] = ...
    optimize.twomey_op(A,b,x_init,1:500,x0);
t.two = toc;

disp('Completed Twomey.');
disp(' ');

eps.two = 100*norm(x0-x_two)/norm(x0);


%% Twomey-Markowski-Buckley
disp('Performing Twomey-Markowski-Buckley...');
tic;
[x_two_mh,Sf_two_mh,out_two_mh] = optimize.twomark_op(...
    A,b,Lb,grid_x,x_init,35,[1e2,1e-5],x0,'Buckley');
t.two_mh = toc;

eps.two_mh = 100*norm(x0-x_two_mh)/norm(x0);



