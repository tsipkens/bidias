
clear;
clc;
close all;


%-- Load colour schemes --------------------------------------------------%
addpath('cmap');
cm_alt = load_cmap('tempo',255);
cm_b = load_cmap('inferno',255);
cm_b = cm_b(40:end,:);
cm_div = load_cmap('RdBu',200);
load('viridis.mat');



%%
%== STEP 1: Generate phantom (x_t) =======================================%
%   High resolution version of the distribution to be projected to coarse
%   grid to generate x.
span_t = [10^-1.5,10^1.5;20,10^3]; % range of mobility and mass

phantom = Phantom('4',span_t);
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
%== STEP 2: Generate A matrix ============================================%
n_b = [14,50]; %[14,50]; %[17,35];
span_b = grid_t.span;
grid_b = Grid(span_b,...
    n_b,'logarithmic'); % grid for data

prop_pma = kernel.prop_pma;
[A_t,sp] = kernel.gen_A_grid(grid_b,grid_t,prop_pma,'Rm',3);
    % generate A matrix based on grid for x_t and b

% sp = kernel.grid2sp(...
%     prop_pma,grid_b,'Rm',3);
% A_alt = kernel.gen_A(sp,grid_b.elements(:,2),...
%     grid_t,prop_pma);

%%
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
%== STEP 3: Generate data ================================================%
b0 = A_t*x_t; % forward evaluate kernel


%-- Corrupt data with noise ----------------------------------------------%
b0(0<1e-10.*max(max(b0))) = 0; % zero very small values of b

Ntot = 1e5;
[b,Lb] = tools.add_noise(b0,Ntot);

figure(5);
colormap(gcf,cm_b);
grid_b.plot2d_marg(b);

figure(20);
grid_b.plot2d_sweep(b,cm_b);



%%
%== STEP 4: Perform inversions ===========================================%
run_inversions_g;
run_inversions_i;



%%
%== STEP 5: Visualize the results ========================================%
% ind = 69; %out_tk1.ind_min;
% x_plot = out_tk1(ind).x;

[~,ind] = min([out_ed_lam.chi]);
% ind = 40;
x_plot = out_ed_lam(ind).x; %out_ed_par(ind).x;


%-- Plot retrieved solution --------------%
figure(10);
colormap(gcf,[cm;1,1,1]);
grid_x.plot2d_marg(x_plot,grid_t,x_t);
% colorbar;
caxis([0,cmax*(1+1/256)]);


%{
figure(13);
grid_x.plot2d_sweep(x_plot,cm);
%}


%-- Plot posterior uncertainties ---------%
% [~,spo] = tools.get_posterior(...
%     A,Lb,out_tk1(ind).lambda.*out_tk1(1).Lpr);

% out = out_ed_par(ind);
% y = [out.lambda,out.ratio,...
%     out.ld,out.corr];
% Gd = [(y(3)/y(2))^2,y(4)*y(3)^2/y(2);y(4)*y(3)^2/y(2),y(3)^2];
%     % y(2) = ratio, y(3) = ld, y(4) = corr
% Lpr = invert.exp_dist_lpr(grid_x.elements(:,2),...
%     grid_x.elements(:,1),out_ed_par(ind).lambda,Gd);

% Lpr = invert.exp_dist_lpr(grid_x.elements(:,2),...
%     grid_x.elements(:,1),out_ed_opt.lambda,out_ed_opt.Gd);
% [~,spo] = tools.get_posterior(...
%     A,Lb,out_ed_opt.lambda.*Lpr);

Lpr = invert.exp_dist_lpr(grid_x.elements(:,2),...
    grid_x.elements(:,1),lambda_ed_lam,Gd);

[~,spo] = tools.get_posterior(...
    A,Lb,lambda_ed_lam.*Lpr);

figure(12);
colormap(gcf,cm_alt);
grid_x.plot2d(spo);
% hold on;
% plot(log10(grid_b.elements(:,2)),...
%     log10(grid_b.elements(:,1)),'.w');
% hold off;
colorbar;
% caxis([0,cmax.*0.3]);


% disp('Computing Bayes factor...');
% [B,F,C] = tools.get_f_tk(A,b,Lb,out_tk0);
% figure(13);
% semilogx([out.lambda],F+C.*grid_b.Ne/grid_x.Ne);
    % extra factor accounts for rank of A


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

%%
%{
%-- Bar plot of results --------------------------------------------------%
figure(30);
chi_names = fieldnames(chi);
chi_vals = zeros(length(chi_names),1);
for ii=1:length(chi_names)
    chi_vals(ii) = chi.(chi_names{ii});
end

bar(chi_vals);
% ylim([0,20]);
% ylim([0,100]);
set(gca,'xticklabel',chi_names);


%%
%{
%-- Bar plot of times ----------------------------------------------------%
figure(40);
t_names = fieldnames(t);
t_vals = zeros(length(t_names),1);
for ii=1:length(t_names)
    t_vals(ii) = mean(t.(t_names{ii}),2);
end

bar(t_vals);
set(gca,'xticklabel',t_names);
set(gca,'yscale','log');



%%
%-- Plot marginal distributions ------------------------------------------%
figure(31);
clf;
dim = 2;

grid_t.plot_marginal(x_t,dim);
grid_x.plot_marginal(...
    {x_Tk1,x_init,x_MART,x_Two,x_TwoMH},dim,x0);



%%
%-- Plot conditional distributions ---------------------------------------%
figure(31);
clf;
dim = 2;
ind_plot = 25;

grid_x.plot_conditional(...
    {x0,x_Tk1,x_init,x_MART,x_Two,x_TwoMH},dim,ind_plot,x0);
%}
%}
