
% MAIN_BAYES  Script associated with Sipkens et al., J. Aerosol Sci. (Submitted).
% Executes and compares Bayesian inversions schemes.
% Author: Timothy Sipkens
%=========================================================================%

clear;
clc;
close all;


%-- Load colour maps -----------------------------------------------------%
addpath('cmap');
cm_alt = load_cmap('BuPu',255);
cm_b = load_cmap('inferno',255);
cm_b = cm_b(40:end,:);
cm_div = load_cmap('RdBu',200);
load('viridis.mat');



%%
%== STEP 1: Generate phantom (x_t) =======================================%
%   High resolution version of the distribution to be projected to coarse
%   grid to generate x.
span_t = [10^-1.5,10^1.5;20,10^3]; % range of mobility and mass

phantom = Phantom('1',span_t);
x_t = phantom.x;
grid_t = phantom.grid;
nmax = max(x_t);
cmax = nmax;

%--  Generate x vector on coarser grid -----------------------------------%
n_x = [50,64]; % number of elements per dimension in x
    % [20,32]; % used for plotting projections of basis functions
    % [40,64]; % used in evaluating previous versions of regularization

grid_x = Grid([grid_t.span],...
    n_x,'logarithmic');
x0 = grid_x.project(grid_t,x_t); % project into basis for x

figure(1);
phantom.plot;
colormap(gcf,[cm;1,1,1]);
caxis([0,cmax*(1+1/256)]);

hold on; % plots mg ridges of phantom
plot(log10(grid_t.edges{2}),...
    log10(phantom.mg_fun(grid_t.edges{2})),'w:');
hold off;



%%
%== STEP 2A: Generate A matrix ===========================================%
%   Note that here a dense kernel is computed for
%   data synthesis in Step 3. 
n_b = [14,50]; %[14,50]; %[17,35];
span_b = grid_t.span;
grid_b = Grid(span_b,...
    n_b,'logarithmic'); % grid for data

prop_pma = kernel.prop_pma;
[A_t,sp] = kernel.gen_grid(grid_b,grid_t,prop_pma,'Rm',3);
    % generate A matrix based on grid for x_t and b

disp('Transform to discretization in <strong>x</strong>...');
B = grid_x.transform(grid_t); % evaluate matrix modifier to transform kernel
A = A_t*B; % equivalent to integration, rebases kernel to grid for x (instead of x_t)
A = sparse(A);
disp('Complete.');
disp(' ');

figure(2);
colormap(gcf,[cm;1,1,1]);
grid_x.plot2d_marg(x0,grid_t,x_t);
caxis([0,cmax*(1+1/256)]);



%%
%== STEP 2b: Generate data ===============================================%
b0 = A_t*x_t; % forward evaluate kernel


%-- Corrupt data with noise ----------------------------------------------%
b0(0<1e-10.*max(max(b0))) = 0; % zero very small values of b

Ntot = 1e5;
[b,Lb] = tools.get_noise(b0,Ntot);

figure(5);
tools.plot2d_scatter(...
    grid_b.elements(:,1),grid_b.elements(:,2),b,cm_b);

figure(6);
colormap(cm_b);
grid_b.plot2d_marg(b);

figure(20);
grid_b.plot2d_sweep(b,cm_b);


[pha_b,Nb] = Phantom.fit2(b,grid_b,2,[0,1.7,0.1,2.3]);

%-- pha_b.Sigma properties --------%
s1b = sqrt(pha_b.Sigma{1}(1,1));
s2b = sqrt(pha_b.Sigma{1}(2,2));
R12b = pha_b.Sigma{1}(1,2)/(s1b*s2b);
Dmb = pha_b.Sigma{1}(1,2)/pha_b.Sigma{1}(2,2); % also s1*R12/s2
%----------------------------------%



%%
%== STEP 3: Perform inversions ===========================================%
run_inversions_g;
run_inversions_i;
run_inversions_j;
run_inversions_k; % time methods



%%
%== STEP 4: Visualize the results ========================================%
ind = out_tk1.ind_min;
x_plot = out_tk1(ind).x;
out = out_tk1;
lambda = out_tk1(ind).lambda;

% [~,ind] = max([out_ed_lam.B]);
% x_plot = out_ed_lam(ind).x;
% out = out_ed_lam;
% lambda = lambda_ed_lam;



%-- Plot retrieved solution --------------%
figure(10);
colormap(gcf,[cm;1,1,1]);
grid_x.plot2d_marg(x_plot,grid_t,x_t);
% colorbar;
caxis([0,cmax*(1+1/256)]);



%-{
%-- Plot posterior uncertainties ------------------------------------%
%-- Tikhonov --%
% [~,spo] = tools.get_posterior(...
%     A,Lb,out_tk1(ind).lambda.*out_tk1(1).Lpr);
% figure(12);
% colormap(gcf,cm_alt);
% grid_x.plot2d(spo);
% colorbar;


%-- Exponential distance --%
%{
if iscell(phantom.Sigma)
    Gd = phantom.Sigma{1};
else
    Gd = phantom.Sigma;
end
Lpr = invert.exp_dist_lpr(Gd,grid_x.elements(:,2),...
    grid_x.elements(:,1));
[~,spo] = tools.get_posterior(...
    A,Lb,lambda.*Lpr);
figure(12);
colormap(gcf,cm_alt);
grid_x.plot2d(spo);
%}



%-- Plot regularization parameter selection schemes ----------------------%
%-{
figure(13);
loglog([out.lambda],[out.chi]); % plot absolute Euclidean error
hold on;
loglog([out.lambda],-([out.B])); % plot Bayes factor
loglog([out.lambda],-([out.F]),'--'); % plot fit
loglog([out.lambda],-([out.C]),'--'); % plot credence
hold off;
%}



figure(10);


