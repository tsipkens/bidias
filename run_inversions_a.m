
% RUN_INVERSIONS_A  Run inversions to find optimal regularization parameters
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

%% 
% Initial guess for iterative schemes
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


%% 
% Least squares
disp('LS inversion ...');
x_length = length(A(1,:));
x_lsq = invert.lsq(Lb*A,Lb*b);
tools.textdone();
disp(' ');
eps.lsq = norm(x0 - x_lsq);


%% 
% Tikhonov (0th) implementation
tools.textheader('Tikhonov (0th) regularization');
tic;
[x_tk0,lambda_tk0,out_tk0] = optimize.tikhonov_op(...
    Lb*A,Lb*b,[1e-2,1e2],0,grid_x,x0,[],'interior-point');
t.tk0 = toc;
disp('Inversion complete.');
disp(' ');

eps.tk0 = norm(x0-x_tk0);


%% 
% Tikhonov (1st) implementation
tools.textheader('Tikhonov (1st) regularization');
tic;
[x_tk1,lambda_tk1,out_tk1] = optimize.tikhonov_op(...
    Lb*A,Lb*b,[1e-2,1e2],1,grid_x,x0,[],'interior-point');
t.tk1 = toc;
disp('Inversion complete.');
disp(' ');

eps.tk1 = norm(x0-x_tk1);


%% 
% Tikhonov (2nd) implementation
tools.textheader('Performing Tikhonov (2nd) regularization');
tic;
[x_tk2,lambda_tk2,out_tk2] = optimize.tikhonov_op(...
    Lb*A,Lb*b,[1e-2,1e2],2,grid_x,x0,[],'interior-point');
t.tk2 = toc;
disp('Inversion complete.');
disp(' ');

eps.tk2 = norm(x0-x_tk2);


%% 
% MART, Maximum entropy regularized solution

tools.textheader('MART');
tic;
[x_mart,iter_mart,out_mart] = ...
    optimize.mart_op(A,b,x_init,1:300,x0);
t.MART = toc;
disp('Inversion complete.');
disp(' ');

eps.mart = norm(x0-x_mart);


%% 
% Twomey
tools.textheader('Twomey');
tic;
[x_two,iter_two,out_two] = ...
    optimize.twomey_op(A,b,x_init,1:500,x0);
t.two = toc;

disp('Completed Twomey.');
disp(' ');

eps.two = norm(x0-x_two);


%% 
% Twomey-Markowski-Buckley
tools.textheader('Twomey-Markowski-Buckley');
tic;
[x_two_mh,Sf_two_mh,out_two_mh] = ...
    optimize.twomark_op(A,b,Lb,grid_x,...
    x_init,35,[1e2,1e-5],x0,'Buckley');
t.two_mh = toc;

eps.two_mh = norm(x0-x_two_mh);



