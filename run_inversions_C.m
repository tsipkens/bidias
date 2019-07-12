

%% Initial guess for iterative schemes
x_init = interp2(grid_b.edges{2}',grid_b.edges{1}',...
reshape(full(b)./(A*ones(size(x0))),grid_b.ne),...
grid_x.elements(:,2),grid_x.elements(:,1));
x_init(isnan(x_init)) = 0;
x_init(isinf(x_init)) = 0;
x_init = sparse(max(0,x_init));
chi.init = norm(x0-x_init);


%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
lambda_Tk0 = 0.1105;
tic;
[x_Tk0,D_Tk0,L_Tk0] = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk0,0,sparse(x0));
t.Tk0 = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk0 = norm(x0-x_Tk0);


%% Tikhonov (1st) implementation
disp('Performing Tikhonov (1st) regularization...');
lambda_Tk1 = 5.3044;
tic;
[x_Tk1,D_Tk1,L_Tk1] = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk1,1,sparse(x0));
t.Tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk1 = norm(x0-x_Tk1);


%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
lambda_Tk2 = 6.0619;
tic;
[x_Tk2,D_Tk2,L_Tk2] = invert.tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk2,2,sparse(x0));
t.Tk2 = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk2 = norm(x0-x_Tk2);


%% Exponential
%{
s1 = 0.5;
s2 = 0.5;
Lex = diag([1/s1,1/s2]);
lambda_exp = 10;

disp('Performing rotated exponential distance regularization...');
tic;
[x_exp,L] = invert.exponential_distance(Lb*A,Lb*b,...
    grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_exp,Lex,x0);
t.exp = toc;
disp('Inversion complete.');
disp(' ');

chi.exp = norm(x0-x_exp);


%% Exponential, rotated

s1 = 1.0;
s2 = 0.1;
dtot = @(d1,d2) sqrt(exp(d1).^2+exp(d2).^2);
theta = -atan2(1,3);%-45/180*pi;%-atan2(3,1);
Lex = diag([1/s1,1/s2])*...
    [cos(theta),-sin(theta);sin(theta),cos(theta)];
lambda_expRot = 10;

disp('Performing rotated exponential distance regularization...');
tic;
[x_expRot,L] = invert.exponential_distance(Lb*A,Lb*b,...
    grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_expRot,Lex,x0);
t.expRot = toc;
disp('Inversion complete.');
disp(' ');

chi.expRot = norm(x0-x_expRot);


%% Exponential, rotated

s1 = 0.8;
s2 = 0.6;
dtot = @(d1,d2) sqrt(exp(d1).^2+exp(d2).^2);
theta = -atan2(1,2);
Lex = diag([1/s1,1/s2])*...
    [cos(theta),-sin(theta);sin(theta),cos(theta)];
lambda_expRot2 = 10;

disp('Performing rotated exponential distance regularization...');
tic;
[x_expRot2,L] = invert.exponential_distance(Lb*A,Lb*b,...
    grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_expRot2,Lex,x0);
t.expRot2 = toc;
disp('Inversion complete.');
disp(' ');

chi.expRot2 = norm(x0-x_expRot2);
%}

