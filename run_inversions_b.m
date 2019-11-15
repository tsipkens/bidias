
% RUN_INVERSIONS_B  Single inversion of each technique using externally defined parameters.
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
x_lsq = invert.lsq(Lb*A,Lb*b);
disp('Inversion complete.');
disp(' ');

disp('Performing LS inversion...');
tic;
x_lsq_nn = invert.lsq(Lb*A,Lb*b,[],'non-neg');
t.lsq = toc;
disp('Inversion complete.');
disp(' ');

% chi.lsq = norm(x0-x_lsq);
chi.lsq = norm(x0-x_lsq_nn);


%% Tikhonov (0th) implementation
% lambda_tk0_hr = tools.perform_hankeraus(out_tk0,A,b,0);
% disp('Performing Tikhonov (0th) regularization (Hanke-Raus)...');
% x_tk0_hr = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk0_hr,0);
% disp('Inversion complete.');
% disp(' ');

disp('Performing Tikhonov (0th) regularization (int.-p.)...');
x_tk0 = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk0,0);
disp('Inversion complete.');
disp(' ');

disp('Performing Tikhonov (0th) regularization (non-neg)...');
x_tk0_nn = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk0,0,[],'non-neg');
disp('Inversion complete.');
disp(' ');

% diff.tk0 = norm(x_tk0-x_tk0_nn);
chi.tk0 = norm(x0-x_tk0_nn);
% chi.tk0_hr = norm(x0-x_tk0_hr);


%% Tikhonov (1st) implementation
% lambda_tk1_hr = tools.perform_hankeraus(out_tk1,A,b,0);
% disp('Performing Tikhonov (1st) regularization (Hanke-Raus)...');
% x_tk1_hr = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk1_hr,1);
% disp('Inversion complete.');
% disp(' ');

disp('Performing Tikhonov (1st) regularization (int.-p.)...');
x_tk1 = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk1,1);
disp('Inversion complete.');
disp(' ');

disp('Performing Tikhonov (1st) regularization (non-neg)...');
x_tk1_nn = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk1,1,[],'non-neg');
disp('Inversion complete.');
disp(' ');

% diff.tk1 = norm(x_tk1-x_tk1_nn);
chi.tk1 = norm(x0-x_tk1_nn);
% chi.tk1_hr = norm(x0-x_tk1_hr);


%% Tikhonov (2nd) implementation
% lambda_tk2_hr = tools.perform_hankeraus(out_tk2,A,b,0);
% disp('Performing Tikhonov (2nd) regularization (Hanke-Raus)...');
% x_tk2_hr = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk2_hr,2);
% disp('Inversion complete.');
% disp(' ');

disp('Performing Tikhonov (2nd) regularization (int.-p.)...');
x_tk2 = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk2,2);
disp('Inversion complete.');
disp(' ');

disp('Performing Tikhonov (2nd) regularization (non-neg)...');
x_tk2_nn = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_tk2,2,[],'non-neg');
disp('Inversion complete.');
disp(' ');

% diff.tk2 = norm(x_tk2-x_tk2_nn);
chi.tk2 = norm(x0-x_tk2_nn);
% chi.tk2_hr = norm(x0-x_tk2_hr);


%% MART, Maximum entropy regularized solution

disp('Performing MART...');
x_MART = invert.mart(A,b,x_init,299);
disp('Inversion complete.');
disp(' ');

chi.mart = norm(x0-x_MART);


%% Twomey
disp('Performing Twomey...');
x_Two = invert.twomey(A,b,x_init,500);
disp('Completed Twomey.');
disp(' ');

chi.two = norm(x0-x_Two);


%% Twomey-Markowski-Buckley
disp('Performing Twomey-Markowski...');
x_two_mh = invert.twomark(A,b,Lb,n_x(1),...
    x_init,35,'Buckley',Sf_two_mh);
disp('Completed Twomey-Markowski.');

chi.two_mh = norm(x0-x_two_mh);



