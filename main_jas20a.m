
% MAIN_JAS20A  Script associated with Sipkens et al., J. Aerosol Sci. (2020a).
% Executes and compares deterministic inversions schemes.
% Author: Timothy Sipkens
%=========================================================================%

clear;
clc;
close all;


%-- Load colour schemes --------------------------------------------------%
addpath cmap;
cm = inferno;
cm = cm(40:end,:);
cm_b = cm;
cm = viridis;



%%
%== STEP 1: Generate phantom (x_t) =======================================%
%   High resolution version of the distribution to be projected to coarse
%   grid to generate x.
span_t = [10^-1.5,10^1.5;10,10^3]; % range of mobility and mass

phantom = Phantom('1',span_t);
x_t = phantom.x;
grid_t = phantom.grid;
nmax = max(x_t);
cmax = nmax;

%== Generate x vector on coarser grid ====================================%
%   This will be used later to gauge accuracy of reconstructions
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
n_b = [14,50]; %[12,50]; %[17,35];
span_b = grid_t.span;
grid_b = Grid(span_b,...
    n_b,'logarithmic'); % should be uniform basis

prop_pma = kernel.prop_pma;
A_t = kernel.gen_grid(grid_b,grid_t,prop_pma,[],'Rm',3);
    % generate A matrix based on grid for x_t and b

disp('Transform to discretization in x...');
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
%== STEP 2b: Generate data ==============================================%
b0 = A_t*x_t; % forward evaluate kernel (high dimension)


%-- Corrupt data with noise ----------------------------------------------%
b0(0<1e-10.*max(max(b0))) = 0; % zero very small values of b

Ntot = 1e5;
[b,Lb] = tools.get_noise(b0,Ntot);


figure(5);
colormap(gcf,cm_b);
grid_b.plot2d_marg(b);

figure(20);
grid_b.plot2d_sweep(b,cm_b);



%%
%== STEP 3: Perform inversions ============================================%
run_inversions_c;

% % optimize the inversion schemes, incurs longer runtimes
% run_inversions_a; % optimize regularization parameter
% run_inversions_b;




%%
%== STEP 4: Plot solution =================================================%
x_plot = x_tk1;

figure(10); % plot reconstruction and marginal distributions
colormap(gcf,[cm;1,1,1]);
grid_x.plot2d_marg(x_plot,grid_t,x_t);
caxis([0,cmax*(1+1/256)]);


%%
%-- Bar plot of Euclidean error ------------------------------------------%
figure(30);
eps_names = fieldnames(eps);
eps_vals = zeros(length(eps_names),1);
for ii=1:length(eps_names)
    eps_vals(ii) = eps.(eps_names{ii});
end

bar(eps_vals);
% ylim([0,20]);
% ylim([0,100]);
set(gca,'xticklabel',eps_names);
