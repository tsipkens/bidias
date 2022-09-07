
% RUN_INVERSIONS_C_PLUS  Shorter implementation with customizable regularization parameters.
%  Adds other methods relative to "run_inversions_c".
%  
%  Author: Timothy Sipkens, 2022-09
%=========================================================================%

tools.textheader('Running inversion suite');
disp(' ');

%% 
% Initial guess for iterative schemes
x_init = interp2(grid_b.edges{2}',grid_b.edges{1}',...
reshape(full(b)./(A*ones(size(x0))),grid_b.ne),...
grid_x.elements(:,2),grid_x.elements(:,1));
x_init(isnan(x_init)) = 0;
x_init(isinf(x_init)) = 0;
x_init = sparse(max(0,x_init));
eps.init = norm(x0-x_init);


%% 
% Least squares
disp('Running LSQ ...');
x_length = length(A(1,:));
tic;
x_lsq = invert.lsq(Lb*A, Lb*b);
t.lsq = toc;
tools.textdone();
disp(' ');
eps.lsq = norm(x0 - x_lsq);


%% 
% Twomey
disp('Running Twomey ...');
tic;
x_two = invert.twomey(A,b,x_init,300,[],[],1);
t.two = toc;
disp(' ');

eps.two = norm(x0-x_two);


%% 
% Twomey-Markowski-Buckley
disp('Running Twomey-Markowski-Buckley ...');
tic;
x_two_mh = invert.twomark(A,b,Lb,n_x(1),...
    x_init,35,'Buckley');
t.twomark = toc;
eps.twomark = norm(x0-x_two_mh);


%% 
% MART, Maximum entropy regularized solution

disp('Running MART ...');
tic;
x_mart = invert.mart(A, b, x_init, 1:300);
t.mart = toc;
tools.textdone();
disp(' ');

eps.mart = norm(x0-x_mart);


%% 
% Tikhonov (0th) implementation
disp('Running Tikhonov (0th) ...');
lambda_tk0 = 0.419941123497942;
tic;
[x_tk0,D_tk0,L_tk0,Gpo_tk0] = invert.tikhonov(...
    Lb*A,Lb*b,lambda_tk0,0,n_x(1),[],'non-neg');
t.tk0 = toc;
tools.textdone();
disp(' ');

eps.tk0 = norm(x0-x_tk0);


%% 
% Tikhonov (1st) implementation
disp('Running Tikhonov (1st) ...');
lambda_tk1 = 0.935436889902617;
tic;
[x_tk1,D_tk1,L_tk1,Gpo_tk1] = invert.tikhonov(...
    Lb*A,Lb*b,lambda_tk1,1,n_x(1),[],'non-neg');
t.tk1 = toc;
tools.textdone();
disp(' ');

eps.tk1 = norm(x0-x_tk1);


%% 
% Tikhonov (2nd) implementation
disp('Running Tikhonov (2nd) ...');
lambda_tk2 = 1.069019204603001;
tic;
[x_tk2,D_tk2,L_tk2] = invert.tikhonov(...
    Lb*A,Lb*b,lambda_tk2,2,n_x(1),[],'non-neg');
t.tk2 = toc;
tools.textdone();
disp(' ');

eps.tk2 = norm(x0-x_tk2);

%% 
% Expectation maximization.
disp('Running expectation-maximization (EM) ...');
tic;
x_em = invert.em(Lb*A, Lb*b, x_init, 3);
t.em = toc;
tools.textdone();
disp(' ');

eps.em = norm(x0-x_em);


%% 
% Exponential distance
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

disp('Running exponential distance ...');
lambda_ed = 1.0826; % found using other run_inversion* scripts
tic;
[x_ed] = ...
    invert.exp_dist(...
    Lb*A,Lb*b,lambda_ed,Gd,...
    grid_x,[]);
t.ed = toc;
tools.textdone();
disp(' ');

eps.ed = norm(x_ed-x0);
