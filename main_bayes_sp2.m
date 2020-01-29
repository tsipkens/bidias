
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
span_t = [10^-1,10^1;...
    10^-1,10^1]; % range of mobility and mass

grid_t = Grid(span_t,...
    [540,540],'logarithmic');
grid_t = grid_t.partial(0,1);
phantom = Phantom('distr-sp2',grid_t);
x_t = phantom.x;
nmax = max(x_t);
cmax = nmax;

%--  Generate x vector on coarser grid -----------------------------------%
n_x = [60,60]; % number of elements per dimension in x

grid_x = Grid([grid_t.span],...
    n_x,'logarithmic');
grid_x = grid_x.partial(0,1);
x0 = grid_x.project(grid_t,x_t); % project into basis for x

phantom_a = Phantom('distr-sp2',grid_x);
x_a = phantom_a.x;

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
n_b = [60,10]; %[14,50]; %[17,35];
span_b = grid_t.span;
grid_b = Grid(span_b,...
    n_b,'logarithmic'); % grid for data
% grid_b = grid_b.partial(0,1);

prop_pma = kernel.prop_pma;
[A_t,sp] = kernel.gen_kernel_grid_c2(grid_b,grid_t,prop_pma,'Rm',3);
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
%== STEP 2b: Generate data ===============================================%
b0 = A_t*x_t; % forward evaluate kernel


%-- Corrupt data with noise ----------------------------------------------%
b0(0<1e-10.*max(max(b0))) = 0; % zero very small values of b

% Note: The total number of particles, Ntot, scales the pdf to equal the 
% number distribution, d2N/dlogA*dlogB. It is also equal to the 
% concentration of particles reduced by the product of the flow 
% rate and overall collection time, that is Ntot = N*Q*t. 
Ntot = 1e5; % converts pdf to counts
[b,Lb] = tools.add_noise(b0,Ntot);

figure(5);
colormap(gcf,cm_b);
grid_b.plot2d(b);

figure(20);
grid_b.plot2d_sweept(b,cm_b);



%%
%== STEP 3: Perform inversions ===========================================%
run_inversions_h;



%%
%== STEP 4: Visualize the results ========================================%
x_plot = x_tk1_a;


%-- Plot retrieved solution --------------%
figure(10);
colormap(gcf,[cm;1,1,1]);
grid_x.plot2d_marg(x_plot,grid_t,x_t);
grid_x.overlay_line([0,0],1,'w:');
% colorbar;
caxis([0,cmax*(1+1/256)]);


%{
figure(13);
grid_x.plot2d_sweep(x_plot,cm);
%}

%{
%-- Plot posterior uncertainties ---------%
%   Tikhonov
[~,spo] = tools.get_posterior(...
    A,Lb,out_tk1(ind).lambda.*out_tk1(1).Lpr);
figure(12);
colormap(gcf,cm_alt);
grid_x.plot2d(spo);
colorbar;


Gd = phantom.Sigma{1};
Lpr = invert.exp_dist_lpr(grid_x.elements(:,2),...
    grid_x.elements(:,1),lambda_ed_lam,Gd);
[~,spo] = tools.get_posterior(...
    A,Lb,lambda_ed_lam.*Lpr);
figure(12);
colormap(gcf,cm_alt);
grid_x.plot2d(spo);


%-- Plot difference to true phantom ------%
figure(11);
diff_plot = (x_plot-x0);%./spo;
scl = max(max(abs(diff_plot)));
grid_x.plot2d(diff_plot);
colormap(cm_div);
caxis([-scl,scl]);


figure(10);


%%
[~,ind_tk1] = min([out_tk1.chi]);
[~,ind_ed] = min([out_ed_lam.chi]);

figure(14);
grid_x.plot2d(log10(out_tk1(ind_tk1).x));
colormap(cm);
colorbar;
caxis([-4,1]);

figure(15);
grid_x.plot2d(log10(out_ed_lam(ind_ed).x));
colormap(cm);
colorbar;
caxis([-4,1]);

figure(16);
grid_x.plot2d(log10(x0));
colormap(cm);
colorbar;
caxis([-4,1]);
%}
