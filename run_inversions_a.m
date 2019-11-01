
% RUN_INVERSIONS_B  Run inversions to find optimal regularization parameters
% Author:           Timothy Sipkens, 2019-05-28
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
chi.init = norm(x0-x_init);
x_init_m = grid_x.marginalize(x_init);


%% Least squares
disp('Performing LS inversion...');
x_length = length(A(1,:));
x_lsq = invert.lsq(Lb*A,Lb*b);
disp('Inversion complete.');
disp(' ');

chi.lsq = norm(x0-x_lsq);


%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
tic;
[x_tk0,lambda_tk0,out_tk0] = optimize.tikhonov(Lb*A,Lb*b,n_x(1),...
    [1e-2,1e2],x0,0,[],'interior-point');
t.tk0 = toc;
disp('Inversion complete.');
disp(' ');

chi.tk0 = norm(x0-x_tk0);


%% Tikhonov (1st) implementation
disp('Performing Tikhonov (1st) regularization...');
tic;
[x_tk1,lambda_tk1,out_tk1] = optimize.tikhonov(Lb*A,Lb*b,n_x(1),...
    [1e-2,1e2],x0,1,[],'interior-point');
t.tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.tk1 = norm(x0-x_tk1);


%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
tic;
[x_tk2,lambda_tk2,out_tk2] = optimize.tikhonov(Lb*A,Lb*b,n_x(1),...
    [1e-2,1e2],x0,2,[],'interior-point');
t.tk2 = toc;
disp('Inversion complete.');
disp(' ');

chi.tk2 = norm(x0-x_tk2);


%% MART, Maximum entropy regularized solution

disp('Performing MART...');
tic;
[x_mart,iter_mart,out_mart] = ...
    optimize.mart(A,b,x_init,1:300,x0);
t.MART = toc;
disp('Inversion complete.');
disp(' ');

chi.mart = norm(x0-x_mart);


%% Twomey
disp('Performing Twomey...');
tic;
[x_two,iter_two,out_two] = ...
    optimize.twomey(A,b,x_init,1:500,x0);
t.two = toc;

disp('Completed Twomey.');
disp(' ');

chi.two = norm(x0-x_two);


%% Twomey-Markowski-Buckley
disp('Performing Twomey-Markowski-Buckley...');
tic;
[x_two_mh,Sf_two_mh,out_two_mh] = ...
    optimize.twomark(A,b,Lb,n_x(1),...
    x_init,35,[1e2,1e-5],x0,'Buckley');
t.two_mh = toc;

chi.two_mh = norm(x0-x_two_mh);



