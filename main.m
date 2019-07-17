
clear;
clc;
close all;


%-- Load colour schemes --------------------------------------------------%
addpath('cmap');
cm = load_cmap('YlGnBu',255);
cm_alt = cm;
load('inferno.mat');
cm = cm(40:end,:);
% load('matter.mat');
cm_b = cm;
load('viridis.mat');


%%
%-- Generate phantom (x_t) -----------------------------------------------%
%   High resolution version of the distribution to be projected to coarse 
%   grid to generate x.

span_t = [10^-1.5,10^1.5;10,10^3]; % range of mobility and mass
    % Hogan lab: -1 -> 1.5

phantom = Phantom('demonstration',span_t);
x_t = phantom.x;
grid_t = phantom.grid;
nmax = max(x_t);
cmax = 5;

figure(1);
phantom.plot;
colormap(gcf,[cm;1,1,1]);
caxis([0,cmax*(1+1/256)]);


%%
%-- Generate A matrix and b vector ---------------------------------------%
n_b = [14,50]; %[12,50]; %[17,35];
span_b = grid_t.span;
grid_b = Grid(span_b,...
    n_b,'logarithmic'); % should be uniform basis

A_t = gen_A(grid_b,grid_t); % generate A matrix based on grid for x_t and b
% load('A_t_v10.mat'); % v10, v11


%%
%--  Generate x vector on coarser grid -----------------------------------%
n_x = [50,64]; % number of elements per dimension in x
    % [20,32]; % used for plotting projections of basis functions
    % [40,64]; % used in evaluating previous versions of regularization

grid_x = Grid([grid_t.span],...
    n_x,'logarithmic');
x0 = grid_x.project(grid_t.edges,x_t); % project into basis for x

disp('Transform to discretization in x...');
B = grid_x.rebase(grid_t); % evaluate matrix modifier to transform kernel
A = A_t*B; % equivalent to integration, rebases kernel to grid for x (instead of x_t)
A = sparse(A);
disp('Complete.');
disp(' ');

figure(2);
colormap(gcf,[cm;1,1,1]);
grid_x.plot2d_marg(x0,grid_t,x_t);
caxis([0,cmax*(1+1/256)]);


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

figure(5);
colormap(gcf,cm_b);
grid_b.plot2d_marg(b);

figure(20);
n2 = floor(grid_b.ne(1));
n3 = floor(length(cm_b(:,1))/n2);
cm_b_mod = cm_b(10:n3:end,:);
set(gca,'ColorOrder',cm_b_mod,'NextPlot','replacechildren');
b_plot_rs = reshape(b,grid_b.ne);
semilogx(grid_b.edges{2},b_plot_rs.*Ntot);
% hold on;
% semilogx(grid_b.edges{2},sum(b_plot_rs).*Ntot,'k');
% hold off;


%%
%-- Perform inversions ---------------------------------------------------%
% run_inversions_A;
% run_inversions_B;
run_inversions_C;


%%
x_plot = x_Tk1;

figure(10);
colormap(gcf,[cm;1,1,1]);
grid_x.plot2d_marg(x_plot,grid_t,x_t);
caxis([0,cmax*(1+1/256)]);

figure(13);
n1 = ceil(grid_x.ne(1)./20);
n2 = floor(grid_x.ne(1)/n1);
n3 = floor(240/n2);
cm_x = cm(10:n3:250,:);
set(gca,'ColorOrder',cm_x,'NextPlot','replacechildren');
x_plot_rs = reshape(x_plot,grid_x.ne);
semilogx(grid_x.edges{2},x_plot_rs(1:n1:end,:));

figure(10);


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

x_t_m = grid_t.marginalize(x_t);
x0_m = grid_x.marginalize(x0);
x_LS_m = grid_x.marginalize(x_LS);
x_Tk0_m = grid_x.marginalize(x_Tk0);
x_Tk1_m = grid_x.marginalize(x_Tk1);
x_init_m = grid_x.marginalize(x_init);
x_MART_m = grid_x.marginalize(x_MART);
x_Two_m = grid_x.marginalize(x_Two);
x_TwoMH_m = grid_x.marginalize(x_TwoMH);

dim = 2;

figure(31);
subplot(3,1,2:3);
semilogx(grid_t.edges{dim},x_t_m{dim},'k');
hold on;
% semilogx(grid_x.edges{dim},x_LS_m{dim});
% semilogx(grid_x.edges{dim},x_Tk0_m{dim});
semilogx(grid_x.edges{dim},x_Tk1_m{dim});
semilogx(grid_x.edges{dim},x_init_m{dim});
semilogx(grid_x.edges{dim},x_MART_m{dim});
semilogx(grid_x.edges{dim},x_Two_m{dim});
semilogx(grid_x.edges{dim},x_TwoMH_m{dim});
hold off;
xlim([min(grid_x.edges{dim}),max(grid_x.edges{dim})]);
% ylim([0,1.4]);

subplot(3,1,1);
semilogx(grid_t.edges{dim},0.*x_t_m{dim},'k');
hold on;
% semilogx(grid_x.edges{dim},x_LS_m{dim}-x0_m{dim});
% semilogx(grid_x.edges{dim},x_Tk0_m{dim}-x0_m{dim});
semilogx(grid_x.edges{dim},x_Tk1_m{dim}-x0_m{dim});
semilogx(grid_x.edges{dim},x_init_m{dim}-x0_m{dim});
semilogx(grid_x.edges{dim},x_MART_m{dim}-x0_m{dim});
semilogx(grid_x.edges{dim},x_Two_m{dim}-x0_m{dim});
semilogx(grid_x.edges{dim},x_TwoMH_m{dim}-x0_m{dim});
hold off;
xlim([min(grid_x.edges{dim}),max(grid_x.edges{dim})]);
ylimits = ylim;
ylim(sign(ylimits).*ceil(abs(ylimits)*20)/20);


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

