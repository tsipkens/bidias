
% MAIN_0  A simple script demonstrating a simple inversion. 


clear;
close all;
clc;

addpath tfer autils cmap;


%== (1) ==================================================================%
span = [20, 800; ...
    10, 600];  % span of grid
ne = [100, 125]; % number of elements in grid for each dimension

% Create an instance of Grid, with logarithmic spacing.
grid_x = Grid(span, ne, 'log');
grid_x.type = {'dm', 'da'};

ut_r = [3,2]; % point in line to cut upper triangle
ut_m = 2; % slope for line to cut upper triangle
lt_r = [1,2]; % point in line to cut lower triangle
lt_m = 2; % slope for line to cut upper triangle
grid_x = grid_x.partial(...
    ut_r, ut_m,...
    lt_r, lt_m); % convert to a partial grid


Gam = corr2cov(log10([1.5, 1.5]'), [1,0.97;0.97,1]);
phantom = Phantom('standard', [], log10([100, 120]), Gam);
x0 = phantom.eval(grid_x);  % evaluate phantom on grid_x

% Plot the phantom.
figure(1);
grid_x.plot2d(x0);
colormap(viridis);  % apply primary colormap



%== (2) ==================================================================%
% Define a new grid for the measurements
span_b = span;
ne_b = [20, 65];
grid_b = Grid(span_b, ne_b, 'log');


% Get properites. 
prop_a = prop_aac();
prop_d = prop_dma();

% Generate the kernel, use default CPMA properties. 
A = kernel.build_grid(grid_b, grid_x, 1:3, ...
    'dma', {prop_d}, 'aac', {prop_a}, 'charger', {});

figure(2);
grid_x.plot2d_marg(A(527,:)); % plot kernel for 527th data point
colormap(inferno);  % apply alternative colormap


%== (3) ==================================================================%
b0 = A * x0; % generate a set of data using the forward model

[b, Lb] = tools.get_noise(b0, 1e5); % corrupt data, assume peak counts ~1e5


%== (4) ==================================================================%
disp('Tikhonov inversion ...');
lambda = 0.2;  % regularization parameter
order = 1;  % order of Tikhonov matrix to be used
x_tk1 = invert.tikhonov(Lb*A, Lb*b, ...
    lambda, order, grid_x);  % Tikhonov solution
tools.textdone();
disp(' ');
disp(' ');
disp(' ');



%== (5) ==================================================================%
figure(4);
grid_x.plot2d(x_tk1);  % plot Tikhonov solution
colormap(viridis);  % apply primary colormap

