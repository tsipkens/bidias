
% MAIN_DMAS  Script for doing tandem DMA inversion.
% This case is for a HTDMA setup also showing growth factors.
% 
% Author: Timothy Sipkens, 2022-06-28

clear;
clc;
close all;


%-- Load colour schemes --------------------------------------------------%
addpath cmap tfer;
cm = inferno;
cm = cm(40:end,:);
cm_b = cm;
cm = viridis;



%%
%== (1) ==================================================================%
%   Phantom and reconstruction grid.
%   High resolution version of the distribution to be projected to coarse
%   grid to generate x.
span_x = [...
    14, 500; ...  % range of mobilities 2
    14, 500];  % range of mobilities 1


%== Generate x vector on coarser grid ====================================%
%   This will be used later to gauge accuracy of reconstructions
n_x = [20,32]; % number of elements per dimension in x
    % [20,32]; % used for plotting projections of basis functions
    % [40,64]; % used in evaluating previous versions of regularization

grid_x = Grid(span_x,...
    n_x,'logarithmic');

Gp = corr2cov(log10([1.4,1.8]), [1,0.9;0.9,1]);
phantom = Phantom('standard', grid_x, ...
    log10([180, 100]), Gp);

figure(1);
phantom.plot;
colormap(viridis);
hold on;
plot(span_x(1,:), span_x(1,:), 'w');
hold off;



%%
%== (2) ==================================================================%
%   Compute kernel.

n_b = [110,115]; %[12,50]; %[17,35];  % size of the data
span_b = grid_x.span;
grid_b = Grid(span_b,...
    n_b,'logarithmic'); % grid for data

prop = prop_dma;

% Generate A matrix based on grid for x and b.
A = kernel.gen_dmas_grid(grid_b, grid_x, prop, prop);
A2 = kernel.gen_dmas_grid(grid_b, grid_x, prop, prop, 0);

idx = round(2 * grid_b.Ne / 3);
grid_b.elements(idx,:)

figure(2);
grid_x.plot2d(grid_x.reshape(A(idx, :)));
colormap(haline);

% Not reneutralized. 
figure(3);
grid_x.plot2d(grid_x.reshape(A2(idx, :)));
colormap(haline);



%%
%== (3) ==================================================================%
%   Generate data.
x0 = phantom.x;

b0 = A * x0;

Ntot = 1e5;
[b,Lb] = tools.get_noise(b0,Ntot);

figure(4);
grid_b.plot2d_marg(b);
colormap(matter);
hold on;
plot(span_x(1,:), span_x(1,:), 'w');
hold off;



%%
%== (4) ==================================================================%
%   Invert.

%-- Tikhonov (1st order) -------------------------------------------------%
disp('Running Tikhonov (1st) ...');
lambda_tk1 = 1.1053;
x_tk1 = invert.tikhonov(...
    Lb*A, Lb*b, lambda_tk1, 1, n_x(1));
tools.textdone();
disp(' ');

eps.tk1_0 = norm(x0 - x_tk1);

figure(5);
grid_x.plot2d_marg(x_tk1);
colormap(viridis);
hold on;
plot(span_x(1,:), span_x(1,:), 'w');
hold off;



%%
%== (5) =======================%
%   Post-process.

% Convert to growth factor as a function of mobility.
[y, grid_y] = tools.mrbc2frac(x0, grid_x, [0.7, 5], 50);

figure(6);
grid_y.plot2d_marg(y);
colormap(tempo);

