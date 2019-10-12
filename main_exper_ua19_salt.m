
clear;
clc;
close all;


%%
%-- Load colour schemes --------------------------------------------------%
addpath('cmap');
cm = load_cmap('YlGnBu',255);
cm_alt = cm;
load('inferno.mat');
cm = cm(40:end,:);
% load('matter.mat');
cm_b = cm;
load('speed.mat');
cm = cm(end:-1:1,:);


%%
%-- Load experimental data -----------------------------------------------%
% load('..\data\Soot-Salt Data UA May 2019\20190509_SMPS\20190509b_SMPS.mat');
    % 0.03% NaCl
    
% load('..\data\Soot-Salt Data UA May 2019\20190509_SMPS\20190509d_SMPS.mat');
    % 0.3% NaCl
    
load('..\data\Soot-Salt Data UA May 2019\20190509_SMPS\20190509f_SMPS.mat');
    % 3% NaCl
    
% load('..\data\Soot-Salt Data UA May 2019\20190509_SMPS\20190509g_SMPS.mat');
    % Dry Propane

%-- Reformat data --------------------------------------------------------%
data = data';
b_max = max(max(data));
b = data(:)./b_max;
edges_b = {data_m,data_d'};
grid_b = Grid(edges_b,...
    [],'logarithmic');
sig = sqrt(data(:))+max(sqrt(data(:))).*0.01; % estimate of noise
Lb = diag(sqrt(1./sig));

figure(5);
colormap(gcf,cm_b);
grid_b.plot2d_marg(b);

figure(20);
n2 = floor(grid_b.ne(1));
n3 = floor(length(cm_b(:,1))/n2);
cm_b_mod = cm_b(10:n3:end,:);
set(gca,'ColorOrder',cm_b_mod,'NextPlot','replacechildren');
b_plot_rs = reshape(b,grid_b.ne);
semilogx(grid_b.edges{2},b_plot_rs.*b_max);
% hold on;
% semilogx(grid_b.edges{2},sum(b_plot_rs).*Ntot,'k');
% hold off;


%%
%-- Generate A and grid_x ------------------------------------------------%
ne_x = [50,64]; % number of elements per dimension in x
    % [20,32]; % used for plotting projections of basis functions
    % [40,64]; % used in evaluating previous versions of regularization

span = [10^-2.5,50;10,10^3];
grid_x = Grid(span,...
    ne_x,'logarithmic');

r_x = grid_x.nodes;
edges_x = grid_x.edges;
n_x = grid_x.ne;

disp('Evaluate kernel...');
A = kernel.gen_A(grid_b,grid_x); % generate A matrix based on grid for x and b


%% 
%-- Perfrom exponential, rotated regularization --------------------------%
s1 = 1.0;
s2 = 0.1;
dtot = @(d1,d2) sqrt(exp(d1).^2+exp(d2).^2);
theta = -atan2(1,2.5);
Lex = diag([1/s1,1/s2])*...
    [cos(theta),-sin(theta);sin(theta),cos(theta)];
lambda_expRot = 1e-3; % 5e-4

disp('Performing rotated exponential distance regularization...');
[x_expRot,L] = invert.exp_dist(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_expRot,Lex);
disp('Inversion complete.');
disp(' ');

x_plot = x_expRot;


%%
%-- Estimate mass-mobility relation --------------------------------------%
figure(40);
colormap(gcf,cm);
grid_x.plot2d(x_plot);
% [Dm,k,rho_100] = grid_x.fit_mass_mob(x_plot,[2.1,-0.1]);
xlabel('log_{10}(d)');
ylabel('log_{10}(m)');

%{
%-- Plots for effective density ------------------------------------------%
[y,grid_rho] = ...
    tools.mass2rho(x_expRot,grid_x);

figure(31);
colormap(gcf,cm_alt);
grid_rho.plot2d_marg(y);
xlabel('log_{10}(d)');
ylabel('log_{10}(\rho_{eff})');

% [m,b] = grid_rho.fit_mass_mob(y,[2,2.8],-0.6);
grid_rho.plot_line_overlay([0,log10(6*k/pi)+9],Dm-3,'w');
rho = 2000; % density of base material
dpe = 10.^((log10(6*k/(pi*rho))+9)/(3-Dm));
%}


