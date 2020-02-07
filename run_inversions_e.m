
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
chi.init = 100*norm(x0-x_init)/norm(x0);
x_init_m = grid_x.marginalize(x_init);


%% Least squares
disp('Performing LS inversion...');
x_length = length(A(1,:));
x_lsq = invert.lsq(Lb*A,Lb*b);
disp('Inversion complete.');
disp(' ');

chi.lsq = 100*norm(x0-x_lsq)/norm(x0);


%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
tic;
[x_tk0,lambda_tk0,out_tk0] = optimize.tikhonov_op(Lb*A,Lb*b,grid_x,...
    [1e-2,1e2],x0,0,[],'interior-point');
t.tk0 = toc;
disp('Inversion complete.');
disp(' ');

chi.tk0 = 100*norm(x0-x_tk0)/norm(x0);


%% Tikhonov (1st) implementation
disp('Performing Tikhonov (1st) regularization...');
tic;
[x_tk1,lambda_tk1,out_tk1] = optimize.tikhonov_op(Lb*A,Lb*b,grid_x,...
    [1e-2,1e2],x0,1,[],'interior-point');
t.tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.tk1 = 100*norm(x0-x_tk1)/norm(x0);


%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
tic;
[x_tk2,lambda_tk2,out_tk2] = optimize.tikhonov_op(Lb*A,Lb*b,grid_x,...
    [1e-2,1e2],x0,2,[],'interior-point');
t.tk2 = toc;
disp('Inversion complete.');
disp(' ');

chi.tk2 = 100*norm(x0-x_tk2)/norm(x0);


%% Tikhonov (0th) implementation L-curve
disp('Performing Tikhonov (0th) regularization (Lcurve method)...');
tic;
order=0;
Lpr0 = tools.Lpr(A,grid_x,order);
Lambdavec=[10^-3,100];

lcurve0=tools.TikhonovLcurve(Lb*A,Lb*b,Lambdavec,Lpr0);
t.tk0_lcurve = toc;
disp('Inversion complete.');
disp(' ');
chi.tk0_lcurve = 100*norm(x0-lcurve0.x)/norm(x0);


%% Tikhonov (1st) implementation L-curve
disp('Performing Tikhonov (0th) regularization (Lcurve method)...');
tic;
order=1;
Lpr0 = tools.Lpr(A,grid_x,order);
Lambdavec=[10^-3,100];

lcurve1=tools.TikhonovLcurve(Lb*A,Lb*b,Lambdavec,Lpr0);
t.tk1_lcurve = toc;
disp('Inversion complete.');
disp(' ');
chi.tk1_lcurve = 100*norm(x0-lcurve1.x)/norm(x0);


%% Tikhonov (2nd) implementation L-curve
disp('Performing Tikhonov (0th) regularization (Lcurve method)...');
tic;
order=2;
Lpr0 = tools.Lpr(A,grid_x,order);
Lambdavec=[10^-3,100];

lcurve2=tools.TikhonovLcurve(Lb*A,Lb*b,Lambdavec,Lpr0);
t.tk2_lcurve = toc;
disp('Inversion complete.');
disp(' ');
chi.tk2_lcurve = 100*norm(x0-lcurve2.x)/norm(x0);


%% MART, Maximum entropy regularized solution

disp('Performing MART...');
tic;
[x_mart,iter_mart,out_mart] = optimize.mart_op(A,b,x_init,1:300,x0);
t.MART = toc;
disp('Inversion complete.');
disp(' ');

chi.mart = 100*norm(x0-x_mart)/norm(x0);


%% Mart with smoothing function
disp('Performing MART with smoothing function...');
tic;
[x_MART_Sm,Sf_MART_Sm,out_MART_Sm] = optimize.MART_Smooth_op(A,b,Lb,grid_x,...
    x_init,35,[1e2,1e-5],x0,'Buckley');
t.MART_Sm = toc;

chi.MART_Sm = 100*norm(x0-x_MART_Sm)/norm(x0);

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

chi.two = 100*norm(x0-x_two)/norm(x0);


%% Twomey-Markowski-Buckley
disp('Performing Twomey-Markowski-Buckley...');
tic;
[x_two_mh,Sf_two_mh,out_two_mh] = ...
    optimize.twomark_op(A,b,Lb,grid_x,...
    x_init,35,[1e2,1e-5],x0,'Buckley');
t.two_mh = toc;

chi.two_mh = 100*norm(x0-x_two_mh)/norm(x0);



