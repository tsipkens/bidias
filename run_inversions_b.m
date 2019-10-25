
% RUN_INVERSIONS_A  Single inversion of each technique using externally defined parameters.
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
tic;
x_lsq = invert.lsq(Lb*A,Lb*b);
t.lsq = toc;
disp('Inversion complete.');
disp(' ');

disp('Performing LS inversion...');
tic;
x_lsq_nn = invert.lsq(Lb*A,Lb*b,[],'non-neg');
t.lsq = toc;
disp('Inversion complete.');
disp(' ');

chi.lsq = norm(x0-x_lsq);
chi.lsq_nn = norm(x0-x_lsq_nn);


%% Tikhonov (0th) implementation
lambda_tk0_hr = tools.perform_hankeraus(out_tk0,A,b,0);
disp('Performing Tikhonov (0th) regularization...');
x_tk0_hr = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk0_hr,0);
disp('Inversion complete.');
disp(' ');

disp('Performing Tikhonov (0th) regularization (int.-p.)...');
x_tk0 = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk0,0);
disp('Inversion complete.');
disp(' ');

disp('Performing Tikhonov (0th) regularization (non-neg)...');
x_tk0_nn = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk0,0,[],'non-neg');
disp('Inversion complete.');
disp(' ');

diff.tk0 = norm(x_tk0-x_tk0_nn);
chi.tk0_hr = norm(x0-x_tk0_hr);
chi.tk0 = norm(x0-x_tk0_nn);


% Tikhonov (1st) implementation
lambda_tk1_hr = tools.perform_hankeraus(out_tk1,A,b,1);
disp('Performing Tikhonov (1st) regularization...');
x_tk1_hr = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk1_hr,1);
disp('Inversion complete.');
disp(' ');

disp('Performing Tikhonov (1st) regularization (int.-p.)...');
x_tk1 = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk1,1);
disp('Inversion complete.');
disp(' ');

disp('Performing Tikhonov (1st) regularization (non-neg)...');
x_tk1_nn = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk1,1,[],'non-neg');
disp('Inversion complete.');
disp(' ');

diff.tk1 = norm(x_tk1-x_tk1_nn);
chi.tk1 = norm(x0-x_tk1_nn);
chi.tk1_hr = norm(x0-x_tk1_hr);


% Tikhonov (2nd) implementation
lambda_tk2_hr = tools.perform_hankeraus(out_tk2,A,b,0);
disp('Performing Tikhonov (2nd) regularization...');
x_tk2_hr = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk2_hr,2);
disp('Inversion complete.');
disp(' ');

disp('Performing Tikhonov (2nd) regularization (int.-p.)...');
x_tk2 = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk2,2);
disp('Inversion complete.');
disp(' ');

disp('Performing Tikhonov (2nd) regularization (non-neg)...');
x_tk2_nn = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk2,2,[],'non-neg');
disp('Inversion complete.');
disp(' ');

diff.tk2 = norm(x_tk2-x_tk2_nn);
chi.tk2 = norm(x0-x_tk2_nn);
chi.tk2_hr = norm(x0-x_tk2_hr);


%% MART, Maximum entropy regularized solution

disp('Performing MART...');
x_MART = invert.mart(A,b,x_init,299);
disp('Inversion complete.');
disp(' ');

chi.MART = norm(x0-x_MART);


%% Twomey
disp('Performing Twomey...');
tic;
x_Two = invert.twomey(A,b,x_init,500);
t.Two = toc;
disp('Completed Twomey.');
disp(' ');

chi.Two = norm(x0-x_Two);


%% Twomey-Markowski-Buckley
disp('Performing Twomey-Markowski...');
tic;
x_TwoMH = invert.twomark(A,b,Lb,n_x(1),...
    x_init,35,'Buckley',1/Sf_TwoMH);
t.TwoMH = toc;
disp('Completed Twomey-Markowski.');

chi.TwoMH = norm(x0-x_TwoMH);



%% Exponential, rotated
%{
s1 = 1.0;
s2 = 0.1;
dtot = @(d1,d2) sqrt(exp(d1).^2+exp(d2).^2);
theta = -atan2(1,3);%-45/180*pi;%-atan2(3,1);
Lex = diag([1/s1,1/s2])*...
    [cos(theta),-sin(theta);sin(theta),cos(theta)];
lambda_expRot = 5;

disp('Performing rotated exponential distance regularization...');
tic;
[x_expRot,L] = invert.exponential_distance(Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_expRot,Lex,x0);
t.expRot = toc;
disp('Inversion complete.');
disp(' ');

chi.expRot = norm(x0-x_expRot);
%}



