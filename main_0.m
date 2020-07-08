
clear;
close all;
clc;

%== STEP 1 ===============================================================%
span = [0.01,100; ...
    10,1000]; % span of grid [min(mass),max(mass); min(mobility),max(mobility)]
ne = [100, 125]; % number of elements in grid for each dimension
grid_x = Grid(span, ne, 'logarithmic'); % create instance of Grid, with logarithmic spacing


ut_r = [2,0.7]; % point in line to cut upper triangle
ut_m = 3; % slope for line to cut upper triangle
lt_r = [2,-0.8]; % point in line to cut lower triangle
lt_m = 3; % slope for line to cut upper triangle
grid_x = grid_x.partial(...
    fliplr(ut_r),ut_m,...
    fliplr(lt_r),lt_m); % convert to a partial grid


phantom = Phantom('4', grid_x); % get Phantom 4 from Sipkens et al. (2020a)
x0 = phantom.x;

figure(1);
grid_x.plot2d_marg(x0);


% define a new grid for the measurements
span_b = span;
ne_b = [20, 65];
grid_b = Grid(span_b, ne_b, 'logarithmic');


%== Step 2A ==============================================================%
prop_pma = kernel.prop_pma % use default CPMA properties (will display in command line)
A = kernel.gen_grid(grid_b, grid_x); % generate the kernel, use default CPMA properties

figure(2);
grid_x.plot2d_marg(A(530,:)); % plot kernel for 530th data point


%== STEP 2b ==============================================================%
b0 = A * x0; % generate a set of data using the forward model

[b, Lb] = tools.get_noise(b0, 1e5); % corrupt data, assume peak counts ~1e5

% plot resultant data as mobility scans at a range of mass setpoint
figure(3);
tools.plot2d_slices(grid_b, b0);
xlabel('log_{10}(d_m)');
ylabel('log_{10}(m_p)');


%== STEP 3 ===============================================================%
disp('Performing Tikhonov regularization...');
lambda = 1; % regularization parameter
order = 1; % order of Tikhonov matrix to be used
x_tk1 = invert.tikhonov(Lb*A, Lb*b, ...
    lambda, order, grid_x); % tikhonov solution
disp('Complete.');
disp(' ');


disp('Performing exponential distance regularization...');
lambda = 1; % regularization parameter
Gd = phantom.Sigma;
x_ed = invert.exp_dist(Lb*A, Lb*b, ...
    lambda, Gd, grid_x); % exponential distance solution
disp('Complete.');
disp(' ');


%== STEP 4 ===============================================================%
figure(4);
subplot(1,2,1);
grid_x.plot2d(x_tk1); % plot Tikhonov solution

subplot(1,2,2);
grid_x.plot2d(x_ed); % plot exponential distance solution

set(gcf, 'Position', [50 300 900 300]); % position and size plot

