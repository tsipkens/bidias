
clear;
close all;
clc;

% load colour schemes
load('cm_inferno.mat');
cm = cm(40:end,:);
cm_inferno = cm;
load('cm_magma.mat');
cm = cm(40:end,:);
cm_magma = cm;
load('cm_plasma.mat');
cm_plasma = cm;
load('cm_viridis.mat');


%%

% load('..\Data\Soot Data FlareNet 18\20180601_E.mat');
% load('..\Data\Soot-Salt Data M9 Flame\data_flameM9_soot+salt_v1.mat');
load('..\Data\Soot-Salt Data UA May 2019\20190509_SMPS\20190509b_SMPS.mat');

data = data';
b_max = max(max(data));
b = data(:)./b_max;
edges_b = {data_m,data_d'};
grid_b = Grid(edges_b,...
    [],'logarithmic');
sig = sqrt(data(:))+max(sqrt(data(:))).*0.01;
Lb = diag(sqrt(1./sig));
%}

figure(3); % plot data, that is b
colormap(gcf,cm_inferno);
grid_b.plot2d(b);

figure(20);
cm_end_shift = 0;
cm_n = length(cm_inferno)-cm_end_shift;
cm_m = (log10(data_m)-log10(0.01))./(log10(20)-log10(0.01));
cm_int = round(cm_n./(length(data_m)+1));
% cm_b = cm_inferno(1:cm_int:cm_n,:);
cm_b = cm_inferno(round(cm_m.*(cm_n-1))+1,:);
% cm_b = cm_inferno(1:22:(end-17),:);
% cm_b = cm_inferno(1:13:(end-5),:);
set(gca,'ColorOrder',cm_b,'NextPlot','replacechildren');
semilogx(grid_b.edges{2},data);
% semilogx(grid_b.edges{2},reshape(b,grid_b.ne));


%%  Generate A and grid_x

ne_x = [50,80]; % number of elements per dimension in x
    % [20,32]; % used for plotting projections of basis functions
    % [40,64]; % used in evaluating previous versions of regularization

span = [10^-2,20;10,10^3];
    % Hogan lab: -1 -> 1.5
grid_x = Grid(span,...
    ne_x,'logarithmic');

r_x = grid_x.nodes;
edges_x = grid_x.edges;
n_x = grid_x.ne;

disp('Evaluate kernel...');
% load('A_3.mat');
A = gen_A(grid_b,grid_x); % generate A matrix based on grid for x and b


%% Exponential, rotated
%-{
s1 = 1.0;
s2 = 0.1;
dtot = @(d1,d2) sqrt(exp(d1).^2+exp(d2).^2);
theta = -atan2(1,2.5);%-45/180*pi;%-atan2(3,1);
Lex = diag([1/s1,1/s2])*...
    [cos(theta),-sin(theta);sin(theta),cos(theta)];
lambda_expRot = 2e-4; % 5e-4

disp('Performing rotated exponential distance regularization...');
tic;
[x_expRot,L] = exponential_distance(Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_expRot,Lex);
t.expRot = toc;
disp('Inversion complete.');
disp(' ');

figure(7);
colormap(gcf,cm);
grid_x.plot2d(log10(x_expRot));
colorbar;
caxis([1,inf]);

figure(11);
marg_dim = 2;
x_expRot_m = grid_x.marginalize(x_expRot);
hold on;
stairs(log10([grid_x.edges{marg_dim},grid_x.span(marg_dim,2)]),[x_expRot_m{marg_dim},0]);
hold off;
%}


%% Tikhonov (1st) implementation

%-{
disp('Performing Tikhonov (1st) regularization...');
lambda_Tk1 = 2e-4;
[x_Tk1,D_Tk1,L_Tk1] = tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk1,1);
disp('Inversion complete.');
disp(' ');

figure(5);
colormap(gcf,cm);
grid_x.plot2d(log10(x_Tk1));
colorbar;
caxis([1,inf]);
%}


%%

%{

%% Least squares
disp('Performing LS inversion...');
A = sparse(A); % no benefit to using full matrix
Lb = sparse(Lb);
x_LS = lsqlin(Lb*A,Lb*b,...
        [],[],[],[],zeros(size(x0)));
disp('Inversion complete.');
Lb = full(Lb);
A = full(A);

figure(4);
% colormap(gcf,[cm;1,1,1]);
grid_x.plot2d(x_LS);
caxis([0,nmax+1e-21]);

chi.LSQ = norm(x0-x_LS)

marg_dim = 2;
figure(11);
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_LS_m = grid_x.marginalize(x_LS);
stairs(log10([grid_x.edges{marg_dim},grid_x.span(marg_dim,2)]),[x_LS_m{marg_dim},0]);
hold off;

%% Tikhonov (0th) implementation
disp('Performing Tikhonov (0th) regularization...');
lambda_Tk0 = 0.1e-3;%0.56;%2.8184e-05;
[x_Tk0,D_Tk0,L_Tk0] = tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk0,0);
% span_Tk = 10.^[2,-2];
% [x_Tk0,lambda_Tk0,out_Tk0] = tikhonov_optimized(Lb*A,Lb*b,n_x(1),span_Tk,x0,0);
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
lambda_Tk1 = 0.1e-3;%2e-3;%5.623;%5.6234e-05;
[x_Tk1,D_Tk1,L_Tk1] = tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk1,1);
% span_Tk = 10.^[2,-2];
% [x_Tk1,lambda_Tk1,out_Tk1] = tikhonov_optimized(Lb*A,Lb*b,n_x(1),span_Tk,x0,1);
disp('Inversion complete.');
disp(' ');

figure(5);
colormap(gcf,cm);
% subplot(1,3,2);
grid_x.plot2d(x_Tk1);
% caxis([0,1200]);

chi.Tk1 = norm(x0-x_Tk1)

marg_dim = 2;
figure(11);
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_Tk1_m = grid_x.marginalize(x_Tk1);
stairs(log10([grid_x.edges{marg_dim},grid_x.span(marg_dim,2)]),[x_Tk1_m{marg_dim},0]);
hold off;

load('..\..\Program - LII\LII Program 3.9\+CENIDE\viridis.mat');
cm_v = cm(10:15:250,:);
set(gca,'ColorOrder',cm_v,'NextPlot','replacechildren');
x_Tk1_rs = reshape(x_Tk1,grid_x.ne);
semilogx(grid_x.edges{2},x_Tk1_rs(1:3:end,:));

load('..\..\Program - LII\LII Program 3.9\+CENIDE\viridis.mat');
cm_v = cm(10:12:250,:);
set(gca,'ColorOrder',cm_v,'NextPlot','replacechildren');
x_Tk1_rs = reshape(x_Tk1,grid_x.ne);
semilogx(grid_x.edges{1},x_Tk1_rs(:,20:2:(end-20))');

%% Tikhonov (2nd) implementation
disp('Performing Tikhonov (2nd) regularization...');
lambda_Tk2 = 1e-3;%4.5;%7.0795e-05;
[x_Tk2,D_Tk2,L_Tk2] = tikhonov(Lb*A,Lb*b,n_x(1),lambda_Tk2,2);
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
% x_MART = mart(Lb*A,Lb*b,x_MART_init,50);
[x_MART,iter_MART,out_MART] = mart_optimized(Lb*A,Lb*b,x_MART_init,1:300,x0);
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
% [x_Two,SIGMA2] = twomey(A_Two,b_Two,x_Two_init,50);
[x_Two,iter_Two,out_Two] = twomey_optimized(A_Two,b_Two,x_Two_init,1:500,x0);

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

% x_TwoMH = twomey_markowski(A_Two,b_Two,Lb,n_x(1),x_Two_init,35,'Buckley',1/10);
[x_TwoMH,Sf_TwoMH,out_TwoMH] = ...
    twomey_markowski_optimized(A_Two,b_Two,Lb,n_x(1),x_Two_init,35,[2,100],x0,'Buckley');

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


%% Statistical approach: uncertainty analysis

t0 = Lb*A;
S_li_inv = t0'*t0;
s_li = sqrt(1./diag(S_li_inv));

figure(50);
colormap(gcf,cm_magma);
% grid_x.plot2d(s_li./max(x0));
% caxis([0,6]);
grid_x.plot2d(log10(s_li./max(x0)));
caxis([-2,log10(6)]);
colorbar;


[~,ind_lambda] = min(out_Tk0.chi);
[~,~,L_Tk0] = tikhonov(Lb*A,Lb*b,n_x(1),out_Tk0.lambda(ind_lambda),1);

S_Tk0_inv = t0'*t0+L_Tk0'*L_Tk0;
s_Tk0 = sqrt(1./diag(S_Tk0_inv));

figure(51);
colormap(gcf,cm_magma);
% grid_x.plot2d(s_Tk0./max(x0));
% caxis([0,6]);
grid_x.plot2d(log10(s_Tk0./max(x0)));
caxis([-2,log10(6)]);
colorbar;

[~,ind_lambda] = min(out_Tk1.chi);
[~,~,L_Tk1] = tikhonov(Lb*A,Lb*b,n_x(1),out_Tk1.lambda(ind_lambda),1);

S_Tk1_inv = t0'*t0+L_Tk1'*L_Tk1;
s_Tk1 = sqrt(1./diag(S_Tk1_inv));

figure(52);
colormap(gcf,cm_magma);
% grid_x.plot2d(s_Tk1./max(x_Tk1));
% caxis([0,0.2]);
grid_x.plot2d(log10(s_Tk1./max(x0)));
caxis([-2,log10(6)]);
colorbar;


[~,ind_lambda] = min(out_Tk2.chi);
[~,~,L_Tk2] = tikhonov(Lb*A,Lb*b,n_x(1),out_Tk2.lambda(ind_lambda),1);

S_Tk2_inv = t0'*t0+L_Tk2'*L_Tk2;
s_Tk2 = sqrt(1./diag(S_Tk2_inv));

figure(53);
colormap(gcf,cm_magma);
% grid_x.plot2d(s_Tk2./max(x_Tk2));
% caxis([0,0.2]);
grid_x.plot2d(log10(s_Tk2./max(x0)));
caxis([-2,log10(6)]);
colorbar;


%% Bar plot of results

figure(30);
chi_names = fieldnames(chi);
chi_vals = zeros(length(chi_names),1);
for ii=1:length(chi_names)
    chi_vals(ii) = chi.(chi_names{ii});
end

bar(chi_vals);
ylim([0,9]);
set(gca,'xticklabel',chi_names);


%% Plot marginal distributions

dim = 2;
figure(31);
% hold on;
semilogx(grid_x.edges{dim},x0_m{dim},'o-');
hold on;
semilogx(grid_x.edges{dim},x_LS_m{dim},'o-');
semilogx(grid_x.edges{dim},x_Tk0_m{dim},'o-');
semilogx(grid_x.edges{dim},x_Tk1_m{dim},'o-');
semilogx(grid_x.edges{dim},x_Tk2_m{dim},'o-');
semilogx(grid_x.edges{dim},x_Two_init_m{dim},'o-');
semilogx(grid_x.edges{dim},x_Two_m{dim},'o-');
semilogx(grid_x.edges{dim},x_TwoMH_m{dim},'o-');
% semilogx(bas_x.edges{dim},x_ME_m{dim},'o-');
semilogx(grid_t.edges{dim},x_t_m{dim},'-k');
hold off;
ylim([0,1.4]);


%% Plot conditional distributions

figure(21);
ind_plot = 25;
x_expRot_rs = reshape(x_expRot,grid_x.ne);
x_Tk1_rs = reshape(x_Tk1,grid_x.ne);
x_TwoMH_rs = reshape(x_TwoMH,grid_x.ne);
x_Two_rs = reshape(x_Two,grid_x.ne);
x_MART_rs = reshape(x_MART,grid_x.ne);
x_Tk0_rs = reshape(x_Tk0,grid_x.ne);
x_LS_rs = reshape(x_LS,grid_x.ne);

%-{
semilogx(grid_x.edges{2},x_LS_rs(ind_plot,:));
hold on;
semilogx(grid_x.edges{2},x_Tk1_rs(ind_plot,:));
semilogx(grid_x.edges{2},x_TwoMH_rs(ind_plot,:));
semilogx(grid_x.edges{2},x_Two_rs(ind_plot,:));
semilogx(grid_x.edges{2},x_MART_rs(ind_plot,:));
% semilogx(grid_x.edges{2},x_Tk0_rs(ind_plot,:));
% semilogx(grid_x.edges{2},x_expRot_rs(ind_plot,:));
hold off;
ylim([0,1500]);
xlim([50,500]);
%}

%{
plot3(1.*ones(length(grid_x.edges{2}),1),log10(grid_x.edges{2}),x_LS_rs(ind_plot,:))
hold on;
plot3(2.*ones(length(grid_x.edges{2}),1),log10(grid_x.edges{2}),x_Tk1_rs(ind_plot,:))
plot3(3.*ones(length(grid_x.edges{2}),1),log10(grid_x.edges{2}),x_TwoMH_rs(ind_plot,:))
plot3(4.*ones(length(grid_x.edges{2}),1),log10(grid_x.edges{2}),x_MART_rs(ind_plot,:))
% plot3(5.*ones(length(grid_x.edges{2}),1),log10(grid_x.edges{2}),x_expRot_rs(ind_plot,:))
hold off;
xlim([0,5]);
zlim([0,1000]);
%}

%}



%% Plots for effective density

[y,rho_vec,dNdlogrho] = ...
    massmob2rhomob(x_expRot,grid_x);

figure(30);
semilogx(rho_vec,dNdlogrho);

figure(31);
% imagesc(log10(grid_x.edges{2}),log10(rho_vec),max(log10(y),max(max(log10(y))).*0.1));
imagesc(log10(grid_x.edges{2}),log10(rho_vec),y);
set(gca,'YDir','normal');
colormap(cm);


