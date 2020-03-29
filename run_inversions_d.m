
% RUN_INVERSIONS_D  Invert multuple times to determine CPU time.
% Author:           Timothy Sipkens, 2019-07-22
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
eps.init = norm(x0-x_init);
x_init_m = grid_x.marginalize(x_init);


for ii=1:20

%% Least squares
disp('Performing LS inversion...');
tic;
x_LSQ = invert.lsq(A,b,'interior-point');
t.LSQ(ii) = toc;
disp('Inversion complete.');
disp(' ');

eps.LSQ = norm(x0-x_LSQ);


%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
tic;
x_tk0 = invert.tikhonov(Lb*A,Lb*b,lambda_tk0,0,n_x(1));
t.tk0(ii) = toc;
disp('Inversion complete.');
disp(' ');

eps.tk0(ii) = norm(x0-x_tk0);


%% Tikhonov (1st) implementation
disp('Performing Tikhonov (1st) regularization...');
tic;
x_tk1 = invert.tikhonov(Lb*A,Lb*b,lambda_tk1,1,n_x(1));
t.tk1(ii) = toc;
disp('Inversion complete.');
disp(' ');

eps.tk1(ii) = norm(x0-x_tk1);


%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
% lambda_tk2 = 8e1;
tic;
x_tk2 = invert.tikhonov(Lb*A,Lb*b,lambda_tk2,2,n_x(1));
t.tk2(ii) = toc;
disp('Inversion complete.');
disp(' ');

eps.tk2(ii) = norm(x0-x_tk2);


%% MART, Maximum entropy regularized solution

disp('Performing MART...');
tic;
x_mart = invert.mart(A,b,x_init,299);
t.mart(ii) = toc;
disp('Inversion complete.');
disp(' ');

eps.mart = norm(x0-x_mart);


%% Twomey
disp('Performing Twomey...');
tic;
x_two = invert.twomey(A,b,x_init,500);
t.two(ii) = toc;
disp('Completed Twomey.');
disp(' ');

eps.two = norm(x0-x_two);


%% Twomey-Markowski-Buckley
disp('Performing Twomey-Markowski...');
tic;
x_two_mh = invert.twomark(A,b,Lb,grid_x,...
    x_init,35,'Buckley',1/Sf_two_mh);
t.two_mh(ii) = toc;
disp('Completed Twomey-Markowski.');

eps.two_mh = norm(x0-x_two_mh);


end



