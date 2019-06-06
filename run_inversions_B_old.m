

%% Least squares
disp('Performing LS inversion...');
tic;
x_LS = lsqlin(Lb*A,Lb*b,...
        [],[],[],[],zeros(size(x0)));
t.LS = toc;
disp('Inversion complete.');

chi.LSQ = norm(x0-x_LS)


%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
% lambda_Tk0 = 1;
tic;
[x_Tk0,D_Tk0,L_Tk0] = tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk0,0);
% span_Tk = 10.^[2,-2];
% [x_Tk0,lambda_Tk0,out_Tk0] = tikhonov_optimized(Lb*A,Lb*b,n_x(1),span_Tk,x0,0);
t.Tk0 = toc;
disp('Inversion complete.');
disp(' ');

figure(5);
colormap(gcf,cm);
% subplot(1,3,1);
grid_x.plot2d(x_Tk0);
% caxis([0,nmax]);

chi.Tk0 = norm(x0-x_Tk0)

marg_dim = 2;
figure(11);
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_Tk0_m = grid_x.marginalize(x_Tk0);
stairs(log10([grid_x.edges{marg_dim},grid_x.span(marg_dim,2)]),[x_Tk0_m{marg_dim},0]);
hold off;


%% Tikhonov (1st) implementation
disp('Performing Tikhonov (1st) regularization...');
% lambda_Tk1 = 17.6351;%2e-3;%5.623;%5.6234e-05;
tic;
[x_Tk1,D_Tk1,L_Tk1] = tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk1,1);
t.Tk1 = toc;
% span_Tk = 10.^[2,-2];
% [x_Tk1,lambda_Tk1,out_Tk1] = tikhonov_optimized(Lb*A,Lb*b,n_x(1),span_Tk,x0,1);
disp('Inversion complete.');
disp(' ');

figure(5);
colormap(gcf,cm);
% subplot(1,3,2);
grid_x.plot2d(x_Tk1);
% caxis([0,nmax]);

chi.Tk1 = norm(x0-x_Tk1)

marg_dim = 2;
figure(11);
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_Tk1_m = grid_x.marginalize(x_Tk1);
stairs(log10([grid_x.edges{marg_dim},grid_x.span(marg_dim,2)]),[x_Tk1_m{marg_dim},0]);
hold off;

load('..\..\Program - LII\LII Program 3.9\+CENIDE\viridis.mat');
cm = cm(10:15:250,:);
set(gca,'ColorOrder',cm,'NextPlot','replacechildren');
x_Tk1_rs = reshape(x_Tk1,grid_x.ne);
semilogx(grid_x.edges{2},x_Tk1_rs(1:3:end,:));

load('..\..\Program - LII\LII Program 3.9\+CENIDE\viridis.mat');
cm = cm(10:12:250,:);
set(gca,'ColorOrder',cm,'NextPlot','replacechildren');
x_Tk1_rs = reshape(x_Tk1,grid_x.ne);
semilogx(grid_x.edges{1},x_Tk1_rs(:,20:2:(end-20))');


%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
% lambda_Tk2 = 1e-3;%4.5;%7.0795e-05;
tic;
[x_Tk2,D_Tk2,L_Tk2] = tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk2,2);
t.Tk2 = toc;
% span_Tk = 10.^[2,-2];
% [x_Tk2,lambda_Tk2,out_Tk2] = tikhonov_optimized(Lb*A,Lb*b,n_x(1),span_Tk,x0,2);
disp('Inversion complete.');
disp(' ');

figure(5);
colormap(gcf,cm);
% subplot(1,3,3);
grid_x.plot2d(x_Tk2);
caxis([0,nmax]);
% contourf(log10(edges_x{2}),log10(edges_x{1}),x_rs_Tk2,100,'LineColor','none');
% imagesc(log10(edges_x{2}),log10(edges_x{1}),x_rs_Tk2);

chi.Tk2 = norm(x0-x_Tk2)

marg_dim = 2;
figure(11);
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_Tk2_m = grid_x.marginalize(x_Tk2);
stairs(log10([grid_x.edges{marg_dim},grid_x.span(marg_dim,2)]),[x_Tk2_m{marg_dim},0]);
hold off;


%% Level set method
%{
x_fun = @(phi) max(phi,0);
x_Tk1_alt = tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk1,1,'algebraic');
fun = @(phi) norm([Lb*A*x_fun(phi)-Lb*b;L_Tk1*phi])^2;
lsq_opt = optimoptions('lsqnonlin','Display','iter','Algorithm','Levenberg-Marquardt');
phi = lsqnonlin(fun,x_Tk1_alt,[],[],lsq_opt);

% fmin_opt = optimset('Display','iter');
% phi = fminsearch(fun,x_Tk1_alt,fmin_opt);

%{
phi0_fun = @(r,lnm,lnd) -sqrt(((r(:,1)-lnm)).^2+((r(:,2)-lnd)).^2)+0.5; % sigmoid distance function (SDF)
phi0 = phi0_fun(log10(bas_x.nodes),0,1.0);

dr = bas_x.dr;
norm_grad = @(phi) sqrt(sum(bas_x.grad(phi).^2,2));
grad2 = @(phi) sqrt(sum(bas_x.grad(norm_grad(phi)).^2,2));
% norm_grad = @(phi) norm(bas_x.grad(phi0));

D1 = Lb*A;
D2 = D1'*D1;
D3 = D1'*Lb*b;
D4 = (L_Tk2')*L_Tk2;

x_fun = @(phi) max(phi,0);
F = @(phi) D2*x_fun(phi)-D3;%+1e-10.*D4*phi;
dphidt = @(phi) F(phi).*norm_grad(phi);
phi0 = x0.*0+phi0;

disp('Performing level set method:');
textbar(0);
phi = phi0;
dt = 1e-9;%max(dr./dphidt(phi0));
Ni = 160;
for ii=1:Ni
    if isnan(sum(dphidt(phi)))
        disp('An error occurred.');
        break
    end
    phi = phi - dt.*dphidt(phi);
    if mod(ii,2)==0; textbar(ii/Ni); end
    
    figure(32);
    plot(dt.*dphidt(phi)+phi);
    pause(0.3);
    
end
disp('Complete.');
disp(' ');
%}

figure(32);
colormap(gcf,cm);
bas_x.plot2d(x_fun(phi));
% bas_x.plot2d(x_fun(phi));
% bas_x.plot2d(int8(x_fun(phi)>0));
hold on;
contour(log10(bas_x.edges{2})+log10(bas_x.edges{2}(end))/2-log10(bas_x.edges{2}(end-1))/2,...
    log10(bas_x.edges{1})+log10(bas_x.edges{1}(end))/2-log10(bas_x.edges{1}(end-1))/2,...
    reshape(x_fun(phi)>0,bas_x.nn),1,'LineColor','White');
% bwb = bwboundaries(reshape(phi,bas_x.nn)>0);
% for kk = 1:numel(bwb)
%     plot(log10(bas_x.edges{2}(bwb{kk}(:,2))),...
%         log10(bas_x.edges{1}(bwb{kk}(:,1))),'r','LineWidth',3);
% end
hold off
%}
%}

%{
%% Exponential
% s1 = 0.1;
% s2 = 0.2;
% lambda_expUnStruc = 1e-9;

s1 = 0.3;
s2 = 0.1;
Lex = diag([1/s1,1/s2]);
lambda_expUnStruc = 1e-5;

disp('Performing exponential distance regularization...');
x_expUnStruc = exponential_distance(Lb*A,Lb*b,bas_x.nodes(:,2),bas_x.nodes(:,1),...
    Lex,lambda_expUnStruc);
disp('Inversion complete.');
disp(' ');

figure(6);
colormap(gcf,cm);
% subplot(1,3,1);
bas_x.plot2d(x_expUnStruc);
caxis([0,nmax]);

chi.expUnStruc = norm(x0-x_expUnStruc)


%% Exponential, structured
s1 = 0.3;
s2 = 0.1;
Lex = diag([1/s1,1/s2]);
sigma = sqrt(max(x_Tk1,0));
sigma = 0.5.*sigma+0.5.*mean(sigma);
lambda_expStruc = 1e-4;

disp('Performing exponential distance (structured) regularization...');
x_expStruc = exponential_distance(Lb*A,Lb*b,bas_x.nodes(:,2),bas_x.nodes(:,1),...
    Lex,lambda_expStruc,[],sigma);
disp('Inversion complete.');
disp(' ');

figure(6);
colormap(gcf,cm);
% subplot(1,3,2);
bas_x.plot2d(x_expStruc);
caxis([0,nmax]);

chi.expStruc = norm(x0-x_expStruc)
%}


%% Exponential, rotated

s1 = 1.0;
s2 = 0.1;
dtot = @(d1,d2) sqrt(exp(d1).^2+exp(d2).^2);
theta = -atan2(1,3);%-45/180*pi;%-atan2(3,1);
Lex = diag([1/s1,1/s2])*...
    [cos(theta),-sin(theta);sin(theta),cos(theta)];
lambda_expRot = 5;

disp('Performing rotated exponential distance regularization...');
tic;
[x_expRot,L] = exponential_distance(Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_expRot,Lex);
t.expRot = toc;
disp('Inversion complete.');
disp(' ');

figure(6);
colormap(gcf,cm);
% subplot(1,3,3);
grid_x.plot2d(x_expRot);
% caxis([0,nmax]);

chi.expRot = norm(x0-x_expRot)

figure(11);
marg_dim = 1;
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_expRot_m = grid_x.marginalize(x_expRot);
stairs(log10([grid_x.edges{marg_dim},grid_x.span(marg_dim,2)]),[x_expRot_m{marg_dim};0]);
hold off;


%% MART, Maximum entropy regularized solution

b_Two = b;
b_Two(b_Two<(1e-5*max(b_Two))) = 0;
b_Two = b_Two;
A_Two = A;

x_MART_init = interp2(grid_b.edges{2}',grid_b.edges{1}',...
    reshape(full(b_Two)./(A_Two*ones(size(x0))),grid_b.ne),...
    grid_x.elements(:,2),grid_x.elements(:,1));
x_MART_init(isnan(x_MART_init)) = 0;
x_MART_init(isinf(x_MART_init)) = 0;
x_MART_init = max(0,x_MART_init);
chi.MART_init = norm(x0-x_MART_init);

disp('Performing MART...');
tic;
x_MART = mart(Lb*A,Lb*b,x_MART_init,300);
% [x_MART,iter_MART,out_MART] = mart_optimized(Lb*A,Lb*b,x_MART_init,1:300,x0);
t.MART = toc;
disp('Inversion complete.');
disp(' ');

figure(9);
colormap(gcf,cm);
grid_x.plot2d(x_MART);
caxis([0,nmax]);

chi.MART = norm(x0-x_MART)

figure(11);
marg_dim = 2;
% plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_MART_m = grid_x.marginalize(x_MART);
stairs(log10([grid_x.edges{marg_dim},grid_x.span(marg_dim,2)]),[x_MART_m{marg_dim},0]);
hold off;


%% Twomey

b_Two = b;
b_Two(b_Two<(1e-5*max(b_Two))) = 0;
b_Two = b_Two;
A_Two = A;

%-- Get initial guess for Twomey -----------------------%
x_Two_init = interp2(grid_b.edges{2}',grid_b.edges{1}',...
    reshape(full(b_Two)./(A_Two*ones(size(x0))),grid_b.ne),...
    grid_x.elements(:,2),grid_x.elements(:,1));
x_Two_init(isnan(x_Two_init)) = 0;
x_Two_init(isinf(x_Two_init)) = 0;
x_Two_init = max(0,x_Two_init);
chi.TwoInit = norm(x0-x_Two_init);
x_Two_init_m = grid_x.marginalize(x_Two_init);

figure(7);
colormap(gcf,cm);
% subplot(1,4,1);
grid_x.plot2d(x_Two_init);
caxis([0,nmax]);

figure(11);
marg_dim = 2;
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_Two_m = grid_x.marginalize(x_Two_init);
stairs(log10([grid_x.edges{marg_dim},grid_x.span(marg_dim,2)]),[x_Two_m{marg_dim},0]);
hold off;

%-- Perform Twomey algorithm ----------------------------%
disp('Performing Twomey...');
tic;
x_Two = twomey(A_Two,b_Two,x_Two_init,500);
% [x_Two,iter_Two,out_Two] = twomey_optimized(A_Two,b_Two,x_Two_init,1:500,x0);
t.Two = toc;

disp('Completed Twomey.');
disp(' ');

x_rs_Two = reshape(x_Two,n_x);

figure(8);
colormap(gcf,cm);
% subplot(1,4,2);
grid_x.plot2d(x_Two);
% caxis([0,nmax]);

chi.Two = norm(x0-x_Two)

figure(11);
marg_dim = 2;
% plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_Two_m = grid_x.marginalize(x_Two);
stairs(log10([grid_x.edges{marg_dim},grid_x.span(marg_dim,2)]),[x_Two_m{marg_dim},0]);
hold off;


%% Twomey-Markowski-Buckley

disp('Performing Twomey-Markowski-Buckley...');

%-- Get initial guess for Twomey ----------------------------%
x_Two_init = interp2(grid_b.edges{2}',grid_b.edges{1}',...
    reshape(full(b_Two)./(A_Two*ones(size(x0))),grid_b.ne),...
    grid_x.elements(:,2),grid_x.elements(:,1));
x_Two_init(isnan(x_Two_init)) = 0;
x_Two_init(isinf(x_Two_init)) = 0;
x_Two_init = max(0,x_Two_init);
chi.TwoInit = norm(x0-x_Two_init);

tic;
x_TwoMH = twomey_markowski(A_Two,b_Two,Lb,n_x(1),x_Two_init,35,'Buckley',0.5);
% [x_TwoMH,Sf_TwoMH,out_TwoMH] = ...
%     twomey_markowski_optimized(A_Two,b_Two,Lb,n_x(1),x_Two_init,35,[2,100],x0,'Buckley');
t.TwoMH = toc;

x_rs_Two_MH = reshape(x_TwoMH,n_x);

figure(10);
clf;
semilogy(out_TwoMH.Sf,out_TwoMH.chi);

figure(8);
colormap(gcf,cm);
% subplot(1,4,4);
grid_x.plot2d(x_TwoMH);
% caxis([0,nmax]);

chi.TwoMH = norm(x0-x_TwoMH)

figure(11);
marg_dim = 2;
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_TwoMH_m = grid_x.marginalize(x_TwoMH);
stairs(log10([grid_x.edges{marg_dim},grid_x.span(marg_dim,2)]),[x_TwoMH_m{marg_dim},0]);
hold off;

