
% GEN_DMAS_GRID  Evaluate kernel/transfer functions for DMA-DMA. 
%  This function exploits the grid structure of the data and 
%  reconstruction domain to speed computation.
% 
%  A = kernel.gen_dmas_grid(GRID_B,GRID_I) uses the data grid, GRID_B,
%  and integration grid, GRID_I, which is generally higher resolution. 
% 
%  A = kernel.gen_dmas_grid(GRID_B,GRID_I,PROP_DMA1) uses the 
%  pre-computed DMA properties in PROP_DMA1 for the first DMA. 
% 
%  A = kernel.gen_dmas_grid(GRID_B,GRID_I,PROP_DMA1,PROP_DMA2) adds an
%  inputs for the properties of the second DMA. 
% 
%  A = kernel.gen_dmas_grid(...,B_NEUT) adds a boolean flag of whether
%  the particles are reneutralized. 
%  
%  ------------------------------------------------------------------------
% 
%  NOTE: Cell arrays are used for Omega_mat and Lambda_mat in order to 
%  allow for the use of sparse matrices, which is necessary to 
%  store information on higher resolutions grids
%  (such as those used for phantoms).
% 
%  AUTHOR: Timothy Sipkens, 2018-11-27

function A = gen_dmas_grid(grid_b, grid_i, prop_dma1, prop_dma2, b_neut, varargin)

if ~exist('prop_dma1','var'); prop_dma1 = []; end
if ~exist('prop_dma2','var'); prop_dma2 = []; end

% By default, consider reneutralization.
if ~exist('b_neut','var'); b_neut = []; end
if isempty(b_neut); b_neut = 1; end

%-- Parse measurement set points (b) -------------------------------------%
r_star = grid_b.elements;  % vector of all setpoint pairs
d1_star = r_star(:,2);  % mobility setpoints no. 1
d2_star = r_star(:,1);  % mobility setpoints no. 2
n_b = grid_b.ne;  % dimensions of data grid
N_b = prod(n_b);  % length of data vector


%-- Generate grid for intergration ---------------------------------------%
n_i = grid_i.ne;  % dimensions of integration grid
N_i = grid_i.Ne;  % length of integration vector

r = grid_i.elements;  % get elements from the grid
d1 = r(:,2);  % mobilities at which to compute the transfer function (not setpoints)
d2 = r(:,1);  % masses at which to compute the transfer function (not setpoints)


%-- Start evaluate kernel ------------------------------------------------%
tools.textheader('Computing DMA-DMA kernel');

%== Evaluate particle charging fractions =================================%
z_vec = (1:3)';  % evaluate charge states 1 -> 3
n_z = length(z_vec);  % length of charge state vector
f_z1 = sparse( ...
    kernel.tfer_charge(d1.*1e-9,z_vec)); % get fraction charged for d1 vector
f_z2 = sparse( ...
    kernel.tfer_charge(d2.*1e-9,z_vec)); % get fraction charged for d2 vector

%== STEP 1: Evaluate DMA transfer function ===============================%
%   Note: The DMA transfer function is 1D (only a function of mobility),
%   which is exploited to speed evaluation. The results is 1 by 3 cell, 
%   with one entry per charge state.
disp(' Computing DMA1 contribution:');
Omega_mat = cell(1,n_z); % pre-allocate for speed, one cell entry per charge state
tools.textbar([0, n_z]);
for kk=1:n_z
    Omega_mat{kk} = kernel.tfer_dma( ...
        grid_b.edges{2} .* 1e-9, ...  % DMA setpoints
        grid_i.edges{2}' .* 1e-9, ...  % points for integration
        z_vec(kk), ...  % integer charge state
        prop_dma1);  % DMA properties
    
    % Remove numerical noise in kernel.
    Omega_mat{kk}(Omega_mat{kk} < (1e-7 .* max(max(Omega_mat{kk})))) = 0;
    
    % Duplicate transfer function for repeated mass in grid_i.
	[~,jj] = max(d1==grid_i.edges{2},[],2);

    if b_neut
        Omega_mat{kk} = f_z1(z_vec(kk),:) .* Omega_mat{kk}(:,jj);
    else
        Omega_mat{kk} = Omega_mat{kk}(:,jj);
    end
    
    tools.textbar([kk, n_z]);
end
disp(' Complete.');
disp(' ');
%=========================================================================%


%== STEP 1: Evaluate DMA transfer function ===============================%
%   Note: The DMA transfer function is 1D (only a function of mobility),
%   which is exploited to speed evaluation. The results is 1 by 3 cell, 
%   with one entry per charge state.
disp(' Computing DMA2 contribution:');
Lambda_mat = cell(1,n_z); % pre-allocate for speed, one cell entry per charge state
tools.textbar([0, n_z]);
for kk=1:n_z
    Lambda_mat{kk} = kernel.tfer_dma( ...
        grid_b.edges{1} .* 1e-9, ...  % DMA setpoints
        grid_i.edges{1}' .* 1e-9, ...  % points for integration
        z_vec(kk), ...  % integer charge state
        prop_dma2);  % DMA properties
    
    % Remove numerical noise in kernel.
    Lambda_mat{kk}(Lambda_mat{kk} < (1e-7 .* max(max(Lambda_mat{kk})))) = 0;
    
    % Duplicate transfer function for repeated mass in grid_i.
	[~,jj] = max(d2==grid_i.edges{1},[],2);

    if b_neut
        Lambda_mat{kk} = (f_z1(z_vec(kk),:) .* Lambda_mat{kk}(:,jj));
    else
        Lambda_mat{kk} = Lambda_mat{kk}(:,jj);
    end

    tools.textbar([kk, n_z]);
end
disp(' Complete.');
disp(' ');
%=========================================================================%


%== SETP 3: Combine to compile kernel ====================================%
disp(' Compiling kernel ...');
K = sparse(N_b,N_i);
for kk=1:n_z
    [~,i1] = max(d2_star==grid_b.edges{1},[],2); % index corresponding to PMA setpoint
    [~,i2] = max(d1_star==grid_b.edges{2},[],2); % index correspondng to DMA setpoint
    
    if b_neut
        K = K + ... % charging contribution
            Lambda_mat{kk}(i1,:) .* ... % PMA contribution
            Omega_mat{kk}(i2,:); % DMA contribution
    else
        K = K + f_z1(z_vec(kk),:) .* ... % charging contribution
            Lambda_mat{kk}(i1,:) .* ... % PMA contribution
            Omega_mat{kk}(i2,:); % DMA contribution
    end
end
tools.textdone();  % print orange DONE
%=========================================================================%

dr_log = grid_i.dr; % area of integral elements in [logm,logd]T space
A = bsxfun(@times,K,dr_log'); % multiply kernel by element area
A = sparse(A); % exploit sparse structure

tools.textheader();


end



