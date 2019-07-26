
% RUN_INVERSIONS_E  Run inversions to find optimal regularization parameters
%                   Works only on Tikhonov, consider posterior covariance
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%


%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
tic;
[x_Tk0,lambda_Tk0,out_Tk0] = invert.optimize_tikhonov(Lb*A,Lb*b,n_x(1),[1e-2,1e2],x0,0);
t.Tk0 = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk0 = norm(x0-x_Tk0);


%% Tikhonov (1st) implementation
disp('Performing Tikhonov (1st) regularization...');
tic;
[x_Tk1,lambda_Tk1,out_Tk1] = invert.optimize_tikhonov(Lb*A,Lb*b,n_x(1),[1e-2,1e2],x0,1);
t.Tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk1 = norm(x0-x_Tk1);


%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
tic;
[x_Tk2,lambda_Tk2,out_Tk2] = invert.optimize_tikhonov(Lb*A,Lb*b,n_x(1),[1e-2,1e2],x0,2);
t.Tk2 = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk2 = norm(x0-x_Tk2);
