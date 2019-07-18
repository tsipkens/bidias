
% RUN_INVERSIONS_B  Run inversions to find optimal regularization parameters
% Author:           Timothy Sipkens, 2019-07-17
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
options = optimoptions('lsqlin',...
    'Algorithm','interior-point','Display','none');
tic;
x_LS = lsqlin(Lb*A,Lb*b,...
    [],[],[],[],zeros(length(x0),1),[],sparse(length(x0),1),options);
t.LS = toc;
disp('Inversion complete.');
disp(' ');

chi.LSQ = norm(x0-x_LS);


%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
% lambda_Tk0 = 1e5;
tic;
[x_Tk0,lambda_Tk0,out_Tk0] = invert.optimize_tikhonov(Lb*A,Lb*b,n_x(1),[1e-2,1e2],x0,0,sparse(length(x0),1));
t.Tk0 = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk0 = norm(x0-x_Tk0);


%% Tikhonov (1st) implementation
disp('Performing Tikhonov (1st) regularization...');
% lambda_Tk1 = 8e1;
tic;
[x_Tk1,lambda_Tk1,out_Tk1] = invert.optimize_tikhonov(Lb*A,Lb*b,n_x(1),[1e-2,1e2],x0,1,sparse(length(x0),1));
t.Tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk1 = norm(x0-x_Tk1);


%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
% lambda_Tk1 = 1e2;
tic;
[x_Tk2,lambda_Tk2,out_Tk2] = invert.optimize_tikhonov(Lb*A,Lb*b,n_x(1),[1e-2,1e2],x0,2,sparse(length(x0),1));
t.Tk2 = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk2 = norm(x0-x_Tk2);


%% MART, Maximum entropy regularized solution

disp('Performing MART...');
tic;
% x_MART = mart(Lb*A,Lb*b,x_init,300);
[x_MART,iter_MART,out_MART] = invert.optimize_mart(A,b,x_init,1:300,x0);
t.MART = toc;
disp('Inversion complete.');
disp(' ');

chi.MART = norm(x0-x_MART);


%% Twomey
%-- Perform Twomey algorithm ----------------------------%
disp('Performing Twomey...');
tic;
% x_Two = twomey(A,b,x_init,500);
[x_Two,iter_Two,out_Two] = invert.optimize_twomey(A,b,x_init,1:500,x0);
t.Two = toc;

disp('Completed Twomey.');
disp(' ');

chi.Two = norm(x0-x_Two);


%% Twomey-Markowski-Buckley

disp('Performing Twomey-Markowski-Buckley...');
tic;
[x_TwoMH,Sf_TwoMH,out_TwoMH] = invert.optimize_twomark(A,b,Lb,n_x(1),...
    x_init,35,[1,1e3],x0,'Buckley');
t.TwoMH = toc;

x_rs_Two_MH = reshape(x_TwoMH,n_x);

chi.TwoMH = norm(x0-x_TwoMH);



