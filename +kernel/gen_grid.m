
% GEN_GRID  Evaluate kernel/transfer functions for DMA-PMA to find A, exploiting grid structure.
% Author: Timothy Sipkens, 2018-11-27
% 
% Notes:
% 	1. This function exploits the grid structure to minimize the number of
%	transfer function evaluations. 
%   2. Cell arrays are used for Omega_mat and Lambda_mat in order to 
%   allow for the use of sparse matrices, which is necessary to 
%   store information on higher resolutions grids
%   (such as those used for phantoms).
%   3. By default, this function uses the analytical PMA transfer function 
%   corresponding to Case 1C from Sipkens et al. (Aerosol Sci. Technol. 2020b).
% 
%-------------------------------------------------------------------------%
% Inputs:
%   grid_b      Grid on which the data exists
%   grid_i      Grid on which to perform integration
%   prop_pma    Structure defining the properties of the PMA
%   varargin    Name-value pairs used in evaluating the PMA tfer. fun.
%=========================================================================%

function [A,sp] = gen_grid(grid_b, grid_i, prop_pma, prop_dma, varargin)

addpath tfer_pma; % add mat-tfer-pma package to MATLAB path
if ~exist('prop_pma','var'); prop_pma = []; end
if isempty(prop_pma); prop_pma = kernel.prop_pma; end
    % import properties of PMA
    % use default properties selected by prop_pma function
if ~exist('prop_dma','var'); prop_dma = []; end

    
%-- Parse measurement set points (b) -------------------------------------%
r_star = grid_b.elements;  % vector of all setpoint pairs
m_star = r_star(:,1);  % mass setpoints
d_star = r_star(:,2);  % mobility setpoints
n_b = grid_b.ne;  % dimensions of data grid
N_b = prod(n_b);  % length of data vector


%-- Generate grid for intergration ---------------------------------------%
n_i = grid_i.ne;  % dimensions of integration grid
N_i = grid_i.Ne;  % length of integration vector

r = grid_i.elements;  % get elements from the grid
m = r(:,1);  % masses at which to compute the transfer function (not setpoints)
d = r(:,2);  % mobilities at which to compute the transfer function (not setpoints)


%-- Start evaluate kernel ------------------------------------------------%
tools.textheader('Computing PMA-DMA kernel');

%== Evaluate particle charging fractions =================================%
z_vec = (1:3)';  % evaluate charge states 1 -> 3
n_z = length(z_vec);  % length of charge state vector
f_z = sparse( ...
    kernel.tfer_charge(d.*1e-9,z_vec)); % get fraction charged for d vector


%== STEP 1: Evaluate DMA transfer function ===============================%
%   Note: The DMA transfer function is 1D (only a function of mobility),
%   which is exploited to speed evaluation. The results is 1 by 3 cell, 
%   with one entry per charge state.
disp('Computing DMA contribution:');
Omega_mat = cell(1,n_z); % pre-allocate for speed, one cell entry per charge state
tools.textbar([0, n_b(2), 0, n_z]);
for kk=1:n_z
    Omega_mat{kk} = sparse(n_b(2),n_i(2));% pre-allocate for speed
    
    for ii=1:n_b(2)  % loop over d_star
        Omega_mat{kk}(ii,:) = kernel.tfer_dma( ...
            grid_b.edges{2}(ii) .* 1e-9, ...  % DMA setpoints
            grid_i.edges{2} .* 1e-9, ...  % points for integration
            z_vec(kk), ...  % integer charge state
            prop_dma);  % DMA properties
        tools.textbar([ii, n_b(2), kk, n_z]);
    end
    
    Omega_mat{kk}(Omega_mat{kk}<(1e-7.*max(max(Omega_mat{kk})))) = 0;
        % remove numerical noise in kernel
        
	[~,jj] = max(d==grid_i.edges{2},[],2);
    Omega_mat{kk} = Omega_mat{kk}(:,jj);
        % duplicate transfer function for repeated mass in grid_i
end
disp('Complete.');
disp(' ');
%=========================================================================%


%== STEP 2: Evaluate PMA transfer function ===============================%
disp('Computing PMA contribution:');

tools.textbar([0, n_b(1), 0, n_z]); % initiate textbar
Lambda_mat = cell(1,n_z); % pre-allocate for speed, one cell entry per charge state
sp = get_setpoint(prop_pma,...  % get PMA setpoints
    'm_star', grid_b.edges{1} .* 1e-18, ...  % mass from the grid
    varargin{:});  % extra name-value pair to specify setpoint

for kk=1:n_z  % loop over the charge state
    Lambda_mat{kk} = sparse(n_b(1), N_i);  % pre-allocate for speed
    
    for ii=1:n_b(1)  % loop over m_star
        
        % Evaluate PMA transfer function.
        Lambda_mat{kk}(ii,:) = kernel.tfer_pma(...
            sp(ii), m.*1e-18, ... 
            d.*1e-9, z_vec(kk), prop_pma)'; 
        
        tools.textbar([ii, n_b(1), kk, n_z]);  % update text progress bar
    end
end
disp(' ');
%=========================================================================%


%== SETP 3: Combine to compile kernel ====================================%
disp('Compiling kernel...');
K = sparse(N_b,N_i);
for kk=1:n_z
    [~,i1] = max(m_star==grid_b.edges{1},[],2); % index corresponding to PMA setpoint
    [~,i2] = max(d_star==grid_b.edges{2},[],2); % index correspondng to DMA setpoint
    
    K = K+f_z(z_vec(kk),:).*... % charging contribution
        Lambda_mat{kk}(i1,:).*... % PMA contribution
        Omega_mat{kk}(i2,:); % DMA contribution
end
%=========================================================================%
disp('Kernel compiled.');

dr_log = grid_i.dr; % area of integral elements in [logm,logd]T space
A = bsxfun(@times,K,dr_log'); % multiply kernel by element area
A = sparse(A); % exploit sparse structure

tools.textheader();


end



