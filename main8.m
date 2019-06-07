
clear;
clc;
close all;

%-- Load colour schemes --------------------------------------------------%
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
%-- Generate phantom (t) -------------------------------------------------%
%-- High resolution version of the distribution to be projected to coarse 
%-- grid to generate x.
main_phantom;

span_t = [10^-1.5,10^1.5;10,10^3];
    % Hogan lab: -1 -> 1.5

[x_t,grid_t,mg] = gen_phantom(phantom,span_t);
nmax = max(x_t);
x_t_m = grid_t.marginalize(x_t);

figure(1); % plot phantom
colormap(gcf,[cm;1,1,1]);
grid_t.plot2d(x_t);
hold on;
plot(log10(grid_t.edges{2}),mg(grid_t.edges{2},1)./log(10));
% plot(log10(grid_t.edges{2}),mg(grid_t.edges{2},2)./log(10));
xlim(log10(grid_t.span(2,:)));
ylim(log10(grid_t.span(1,:)));
hold off;


%-- Generate A matrix and b vector ---------------------------------------%
n_b = [14,50];%[12,50]; %[17,35];
span_b = grid_t.span;%[bas_x.edges{1}(1),bas_x.edges{1}(end);bas_x.edges{2}(1),bas_x.edges{2}(end)];
grid_b = Grid(span_b,...
    n_b,'logarithmic'); % should be uniform basis

A_t = gen_A(grid_b,grid_t); % generate A matrix based on grid for x_t and b
% load('A_t_v10.mat'); % v10, v11


%%
%--  Generate x vector on coarser grid -----------------------------------%
ne_x = [50,64]; % number of elements per dimension in x
    % [20,32]; % used for plotting projections of basis functions
    % [40,64]; % used in evaluating previous versions of regularization

grid_x = Grid([grid_t.span],...
    ne_x,'logarithmic');
x0 = grid_x.project(grid_t.edges,x_t); % project into basis for x
x0 = x0(:);

r_x = grid_x.nodes;
edges_x = grid_x.edges;
n_x = grid_x.ne;

disp('Transform to discretization in x...');
K0 = grid_x.phi(grid_t); % transform based on chosen basis functions in x
A = A_t*K0; % equivalent to integration
A = sparse(A);
disp('Complete.');
disp(' ');

%{
figure(2); % plot projected phantom, that is x
colormap(gcf,cm);
grid_x.plot2d(x0);
caxis([0,nmax]);

figure(3);
marg_dim = 1;
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x0_m = grid_x.marginalize(x0);
plot(log10(grid_x.edges{marg_dim}),x0_m{marg_dim},'-o');
hold off;

figure(4);
marg_dim = 2;
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x0_m = grid_x.marginalize(x0);
plot(log10(grid_x.edges{marg_dim}),x0_m{marg_dim},'-o');
hold off;
%}


%%
%-- Generate data --------------------------------------------------------%
b0 = A_t*x_t; % forward evaluate kernel


%-- Corrupt b with noise -------------------------------------------------%
b0(0<1e-10.*max(max(b0))) = 0; % zero very small values of b

Ntot = 1e5;
theta = 1/Ntot;
gamma = max(sqrt(theta.*b0)).*1e-8; % underlying Gaussian noise
% 13 = 1e5, 1e-8; 13b = 1e4, 1e-5;
Sigma = sqrt(theta.*b0+gamma^2); % sum up Poisson and Gaussian noise
Lb = sparse(1:grid_b.Ne,1:grid_b.Ne,1./Sigma,grid_b.Ne,grid_b.Ne);
rng(0);
epsilon = Sigma.*randn(size(b0));
b = sparse(b0+epsilon); % add noise
b = max(b,0); % remove negative values

figure(5); % plot data, that is b
colormap(gcf,cm_inferno);
grid_b.plot2d(b);


figure(20);
n2 = floor(grid_b.ne(1));
n3 = floor(length(cm_inferno(:,1))/n2);
cm_b = cm_inferno(10:n3:end,:);
set(gca,'ColorOrder',cm_b,'NextPlot','replacechildren');
b_plot_rs = reshape(b,grid_b.ne);
semilogx(grid_b.edges{2},b_plot_rs.*Ntot);


% figure(21);
% set(gca,'ColorOrder',cm_b,'NextPlot','replacechildren');
% plot3(grid_b.edges{2},...
%     ones(grid_b.ne(2),1)*(grid_b.edges{1}),...
%     b_plot_rs.*Ntot);
% set(gcf,'Renderer','zbuffer');
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% view([4,20]);


marg_dim = 2;
figure(6);
b_m = grid_b.marginalize(b);
plot(log10(grid_b.edges{marg_dim}),b_m{marg_dim},'o-k');


%%
%-- Perform inversions ---------------------------------------------------%
run_inversions_A;
% run_inversions_B;



%%
x_plot = x_Tk1;

figure(10);
colormap(gcf,[cm;1,1,1]);
grid_x.plot2d(x_plot);
caxis([0,1*(1+1/256)]);
% caxis([0,0.85*(1+1/256)]);
% caxis([0,1.25*(1+1/256)]);

figure(11);
marg_dim = 1;
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_plot_m = grid_x.marginalize(x_plot);
% plot(log10(grid_x.edges{marg_dim}),x_plot_m{marg_dim},'.-');
stairs(log10(grid_x.nodes{marg_dim}),...
    [x_plot_m{marg_dim};0]);
% stairs(log10(grid_x.nodes{marg_dim}),...
%     [x_temp{marg_dim};0]);
hold off;

figure(12);
marg_dim = 2;
plot(log10(grid_t.edges{marg_dim}),x_t_m{marg_dim},'k');
hold on;
x_plot_m = grid_x.marginalize(x_plot);
% plot(log10(grid_x.edges{marg_dim}),x_plot_m{marg_dim});
stairs(log10(grid_x.nodes{marg_dim}),...
    [x_plot_m{marg_dim},0]);
% stairs(log10(grid_x.nodes{marg_dim}),...
%     [x_temp{marg_dim},0]);
hold off;

figure(13);
load('cm_viridis.mat');
n1 = ceil(grid_x.ne(1)./20);
n2 = floor(grid_x.ne(1)/n1);
n3 = floor(240/n2);
cm_b = cm(10:n3:250,:);
set(gca,'ColorOrder',cm_b,'NextPlot','replacechildren');
x_plot_rs = reshape(x_plot,grid_x.ne);
semilogx(grid_x.edges{2},x_plot_rs(1:n1:end,:));

% x_temp = x_plot_m;

figure(10);

%% Statistical approach: uncertainty analysis
%{
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
%}

%% Bar plot of results

figure(30);
chi_names = fieldnames(chi);
chi_vals = zeros(length(chi_names),1);
for ii=1:length(chi_names)
    chi_vals(ii) = chi.(chi_names{ii});
end

bar(chi_vals);
ylim([0,4]);
set(gca,'xticklabel',chi_names);

%% Bar plot of times

figure(40);
t_names = fieldnames(t);
t_vals = zeros(length(t_names),1);
for ii=1:length(t_names)
    t_vals(ii) = mean(t.(t_names{ii}),2);
end

bar(t_vals);
set(gca,'xticklabel',t_names);
set(gca,'yscale','log');


%% Plot marginal distributions

x_LS_m = grid_x.marginalize(x_LS);
x_Tk0_m = grid_x.marginalize(x_Tk0);
x_Tk1_m = grid_x.marginalize(x_Tk1);
x_init_m = grid_x.marginalize(x_init);
x_MART_m = grid_x.marginalize(x_MART);
x_Two_m = grid_x.marginalize(x_Two);
x_TwoMH_m = grid_x.marginalize(x_TwoMH);

dim = 2;
figure(31);
semilogx(grid_t.edges{dim},x_t_m{dim},'k');
hold on;
% semilogx(grid_x.edges{dim},x_LS_m{dim});
semilogx(grid_x.edges{dim},x_Tk0_m{dim});
semilogx(grid_x.edges{dim},x_Tk1_m{dim});
semilogx(grid_x.edges{dim},x_init_m{dim});
semilogx(grid_x.edges{dim},x_MART_m{dim});
semilogx(grid_x.edges{dim},x_Two_m{dim});
semilogx(grid_x.edges{dim},x_TwoMH_m{dim});
hold off;
% ylim([0,1.4]);


%% Plot conditional distributions

figure(21);
ind_plot = 25;
x_expRot_rs = reshape(x_expRot,grid_x.ne);
x_plot_rs = reshape(x_Tk1,grid_x.ne);
x_TwoMH_rs = reshape(x_TwoMH,grid_x.ne);
x_Two_rs = reshape(x_Two,grid_x.ne);
x_MART_rs = reshape(x_MART,grid_x.ne);
x_Tk0_rs = reshape(x_Tk0,grid_x.ne);
x_LS_rs = reshape(x_LS,grid_x.ne);

%-{
semilogx(grid_x.edges{2},x_LS_rs(ind_plot,:));
hold on;
semilogx(grid_x.edges{2},x_plot_rs(ind_plot,:));
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

