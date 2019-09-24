
% RUN_INVERSIONS_F  Runs exponential rotated regularization (orig. main_bayes)
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

% Exponential rotated
s1 = 2.0;
s2 = 0.4;
dtot = @(d1,d2) sqrt(exp(d1).^2+exp(d2).^2);
theta = -atan2(1,2.5);
Lex = diag([1/s1,1/s2])*...
    [cos(theta),-sin(theta);sin(theta),cos(theta)];
lambda_expRot = 1.8; % 5e-4

disp('Performing rotated exponential distance regularization...');
[x_expRot,L] = invert.exp_dist(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_expRot,Lex);
disp('Inversion complete.');
disp(' ');

