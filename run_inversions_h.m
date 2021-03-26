
% RUN_INVERSIONS_H  Run exponential distance and 1st-order Tikhonov regularization.
% Author: Timothy Sipkens, 2020-04-06
%=========================================================================%


%-- Tikhonov (1st order) -------------------------------------------------%
tools.textheader('Tikhonov (1st) regularization');
lambda_tk1 = 1.1053; % found using other run_inversion* scripts
x_tk1 = invert.tikhonov(...
    Lb*A,Lb*b,lambda_tk1,1,n_x(1));
disp('Inversion complete.');
disp(' ');

eps.tk1_0 = norm(x0-x_tk1);




%-- Exponential distance approach ----------------------------------------%
Gd = phantom.Sigma(:,:,1);
if isempty(Gd) % for Phantom 3
    [~,Gd] = phantom.p2cov(phantom.p(2),phantom.modes(2));
end

%-- Gd properties -----------------%
l1 = sqrt(Gd(1,1));
l2 = sqrt(Gd(2,2));
R12 = Gd(1,2)/(l1*l2);
Dm = Gd(1,2)/Gd(2,2); % s1*R12/s2
%----------------------------------%

tools.textheader('Exponential distance regularization');
lambda_ed = 1.0826; % found using other run_inversion* scripts
[x_ed] = ...
    invert.exp_dist(...
    Lb*A,Lb*b,lambda_ed,Gd,...
    grid_x,[]);
disp('Inversion complete.');
disp(' ');

eps.ed_0 = norm(x_ed-x0);


