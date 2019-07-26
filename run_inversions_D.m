
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
chi.init = norm(x0-x_init);
x_init_m = grid_x.marginalize(x_init);


for ii=1:20

%% Least squares
disp('Performing LS inversion...');
tic;
x_LSQ = invert.lsq(A,b,'interior-point');
t.LSQ(ii) = toc;
disp('Inversion complete.');
disp(' ');

chi.LSQ = norm(x0-x_LSQ);


%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
tic;
x_Tk0 = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk0,0);
t.Tk0(ii) = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk0(ii) = norm(x0-x_Tk0);


%% Tikhonov (1st) implementation
disp('Performing Tikhonov (1st) regularization...');
tic;
x_Tk1 = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk1,1);
t.Tk1(ii) = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk1(ii) = norm(x0-x_Tk1);


%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
% lambda_Tk2 = 8e1;
tic;
x_Tk2 = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk2,2);
t.Tk2(ii) = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk2(ii) = norm(x0-x_Tk2);


%% MART, Maximum entropy regularized solution

disp('Performing MART...');
tic;
x_MART = invert.mart(A,b,x_init,299);
t.MART(ii) = toc;
disp('Inversion complete.');
disp(' ');

chi.MART = norm(x0-x_MART);


%% Twomey
disp('Performing Twomey...');
tic;
x_Two = invert.twomey(A,b,x_init,500);
t.Two(ii) = toc;
disp('Completed Twomey.');
disp(' ');

chi.Two = norm(x0-x_Two);


%% Twomey-Markowski-Buckley
disp('Performing Twomey-Markowski...');
tic;
x_TwoMH = invert.twomark(A,b,Lb,n_x(1),...
    x_init,35,'Buckley',1/Sf_TwoMH);
t.TwoMH(ii) = toc;
disp('Completed Twomey-Markowski.');

chi.TwoMH = norm(x0-x_TwoMH);


end



