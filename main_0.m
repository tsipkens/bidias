
% MAIN_0  A simple script demonstrating a simple inversion. 


clear;
close all;
clc;

% Get colormap.
if length(dir('cmap')) == 2  % if `cmap` not downloaded
    cm = parula;  % use Matlab default
    cma = gray;  % use grey as alternate colormap
else
    addpath cmap;
    cm = viridis;  % otherwise use viridis
    cma = inferno;  % use inferno as alternate colormap
    cma = cma(30:end, :);  % adjust dark values
end



%== (1) ==================================================================%
span = [0.01, 100; ...
    10, 1000];  % span of grid
ne = [100, 125]; % number of elements in grid for each dimension

% Create an instance of Grid, with logarithmic spacing.
grid_x = Grid(span, ne, 'log');


ut_r = [0.5,2]; % point in line to cut upper triangle
ut_m = 3; % slope for line to cut upper triangle
lt_r = [-1.2,2]; % point in line to cut lower triangle
lt_m = 3; % slope for line to cut upper triangle
grid_x = grid_x.partial(...
    ut_r, ut_m,...
    lt_r, lt_m); % convert to a partial grid


phantom = Phantom('4');  % get Phantom 4 from Sipkens et al. (2020a)
x0 = phantom.eval(grid_x);  % evaluate phantom on grid_x

% Plot the phantom.
figure(1);
grid_x.plot2d(x0);
colormap(cm);  % apply primary colormap



%== (2) ==================================================================%
% Define a new grid for the measurements
span_b = span;
ne_b = [20, 65];
grid_b = Grid(span_b, ne_b, 'log');


% Use default CPMA and DMA properties (will display in command line).
addpath('tfer');
prop_p = prop_pma()
prop_d = prop_dma()

% Generate the kernel, use default CPMA properties. 
A = kernel.build_grid(grid_b, grid_x, 1:3, ...
    'pma', {prop_p}, 'dma', {prop_d}, 'charger', {});

figure(2);
grid_x.plot2d_marg(A(527,:)); % plot kernel for 527th data point
colormap(cma);  % apply alternative colormap



%== (3) ==================================================================%
b0 = A * x0; % generate a set of data using the forward model

[b, Lb] = tools.get_noise(b0, 1e5); % corrupt data, assume peak counts ~1e5

% plot resultant data as mobility scans at a range of mass setpoint
figure(3);
opts.f_lines = 1;
tools.plot2d_patch(grid_b, b0, cma, [], opts);
xlabel('log_{10}(d_m)');
ylabel('log_{10}(m_p)');



%== (4) ==================================================================%
disp('Tikhonov inversion ...');
lambda = 1;  % regularization parameter
order = 1;  % order of Tikhonov matrix to be used
x_tk1 = invert.tikhonov(Lb*A, Lb*b, ...
    lambda, order, grid_x);  % Tikhonov solution
tools.textdone();
disp(' ');
disp(' ');


disp('Exponential distance inversion ...');
lambda = 1; % regularization parameter
Gd = phantom.Sigma;
x_ed = invert.exp_dist( ...
    Lb*A, Lb*b, ...
    lambda, Gd, grid_x);  % exponential distance solution
tools.textdone();
disp(' ');
disp(' ');



%== (5) ==================================================================%
figure(4);
subplot(1,2,1);
grid_x.plot2d(x_tk1);  % plot Tikhonov solution
colormap(cm);  % apply primary colormap

subplot(1,2,2);
grid_x.plot2d(x_ed);  % plot exponential distance solution
colormap(cm);

set(gcf, 'Position', [50 300 900 300]); % position and size plot

