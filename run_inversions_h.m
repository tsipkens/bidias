
% RUN_INVERSIONS_H  Re-run select exponential rotated and Tikhonov regularization
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

%-- Tikhonov (1st order) -----%
disp('Performing Tikhonov (1st) regularization...');
tic;
x_tk1_nn = invert.tikhonov(...
    Lb*A,Lb*b,n_x(1),lambda_tk1,1,[],'non-neg');
t.tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.tk1_nn = norm(x0-x_tk1_nn);


%-- Exponential, uncorrelated ------%
lambda_exp = 1.1;
s1 = 0.4;
s2 = 0.2;
Gd = diag([s1^2,s2^2]);

disp('Performing exponential distance regularization...');
x_exp_a = invert.exp_dist(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_exp,Gd);
disp('Inversion complete.');
disp(' ');

chi.exp_a = norm(x0-x_exp_a);


%-- Exponential rotated ------%
lambda_exp_rot = 1.4;
s1 = 0.8;
s2 = 0.3;
Dm = 1.92; % used to specify correlation
R12 = Dm*s2/s1; % correlation
Gd = [s1^2,R12*(s1*s2);R12*(s1*s2),s2^2];

[V,~] = eig(Gd);
the = atan(V(1,2)/V(2,2))/pi*180;

disp('Performing rotated exponential distance regularization...');
x_exp_rot_a = invert.exp_dist(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_exp_rot,Gd2);
disp('Inversion complete.');
disp(' ');

chi.exp_rot_a = norm(x0-x_exp_rot_a);



