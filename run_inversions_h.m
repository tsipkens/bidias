
% RUN_INVERSIONS_H  Re-run select exponential rotated and Tikhonov regularization
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

%-- Tikhonov (1st order) -----%
lambda_tk1 = 1.4;

disp('Performing Tikhonov (1st) regularization...');
tic;
x_tk1_a = invert.tikhonov(...
    Lb*A,Lb*b,grid_x,lambda_tk1,1,[],'non-neg');
t.tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.tk1_nn = norm(x0-x_tk1_a);


%-- Exponential rotated ------%
lambda_exp_rot = 1.4;
s1 = 0.4;
s2 = s1;
R12 = 0.96; % correlation
Gd2 = [s1^2,R12*(s1*s2);R12*(s1*s2),s2^2];

[V,~] = eig(Gd2);
the = atan(V(1,2)/V(2,2))/pi*180;

disp('Performing exponential distance regularization...');
x_exp_rot_a = invert.exp_dist(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_exp_rot,Gd2);
disp('Inversion complete.');
disp(' ');

chi.exp_a = norm(x0-x_exp_rot_a);


