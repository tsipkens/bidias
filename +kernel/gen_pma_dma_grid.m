
% GEN_PMA_DMA_GRID  Evaluate kernel/transfer functions for PMA-DMA. 
%  Relative to kernel.gen_pma_dma, this function exploits the grid 
%  structure of the data and reconstruction domain to speed computation.
% 
%  A = kernel.gen_pma_dma_grid(GRID_B,GRID_I) uses the data grid, GRID_B,
%  and integration grid, GRID_I, which is generally higher resolution. 
% 
%  A = kernel.gen_pma_dma_grid(GRID_B,GRID_I,PROP_PMA) uses the 
%  pre-computed particle mass analyzer properties specified by PROP_PMA. 
%  This structure can be generated manually or with the help of the 
%  kernel.prop_pma(...) function. 
% 
%  A = kernel.gen_pma_dma_grid(GRID_B,GRID_I,PROP_PMA,PROP_DMA) uses the 
%  pre-computed differential mobility analyzer properties specified by 
%  PROP_DMA. This structure can be generated manually or with the help of 
%  the kernel.prop_dma(...) function. PROP_PMA can be excluded by supplying
%  an empty input, i.e., PROP_PMA = []. 
% 
%  A = kernel.gen_pma_dma_grid(...,'name',value) specifies a name-value
%  pairs that is passed on as input to the get_setpoint(...) function
%  from the tfer_pma folder. This specifies one quantity other than the
%  setpoint mass that is used to constrain the PMA operating point. For
%  example, {'Rm',10} specifies a resolution PMA resolution of 10, which is
%  used for all of the setpoints. 
%  If not provided, the function uses the default in the get_setpoint(...)
%  function in the tfer_pma folder.
% 
%  [A,SP] = kernel.gen_pma_dma_grid(...) also outputs the PMA setpoint
%  structure (associated with the mat_tfer_pma submodule) for the given
%  grid.
%  
%  ------------------------------------------------------------------------
% 
%  NOTE: Cell arrays are used for Omega_mat and Lambda_mat in order to 
%  allow for the use of sparse matrices, which is necessary to 
%  store information on higher resolutions grids
%  (such as those used for phantoms).
% 
%  AUTHOR: Timothy Sipkens, 2018-11-27

function [A, sp] = gen_pma_dma_grid(grid_b, grid_i, prop_pma, prop_dma, varargin)

% Add mat-tfer-pma package to MATLAB path.
addpath tfer_pma;

% If not given, import default properties of PMA, 
% as selected by prop_pma function.
if ~exist('prop_pma','var'); prop_pma = []; end
if isempty(prop_pma); prop_pma = kernel.prop_pma; end

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
    Omega_mat{kk} = sparse(n_b(2),n_i(2));  % pre-allocate for speed
    
    for ii=1:n_b(2)  % loop over d_star
        Omega_mat{kk}(ii,:) = kernel.tfer_dma( ...
            grid_b.edges{2}(ii) .* 1e-9, ...  % DMA setpoints
            grid_i.edges{2} .* 1e-9, ...  % points for integration
            z_vec(kk), ...  % integer charge state
            prop_dma);  % DMA properties
        tools.textbar([ii, n_b(2), kk, n_z]);
    end
    
    % Remove numerical noise in kernel.
    Omega_mat{kk}(Omega_mat{kk} < (1e-7 .* max(max(Omega_mat{kk})))) = 0;
    
    % Duplicate transfer function for repeated mass in grid_i.
	[~,jj] = max(d==grid_i.edges{2},[],2);
    Omega_mat{kk} = Omega_mat{kk}(:,jj);
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
        Lambda_mat{kk}(ii,:) = kernel.tfer_pma( ...
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
disp('Kernel compiled.');
%=========================================================================%

dr_log = grid_i.dr; % area of integral elements in [logm,logd]T space
A = bsxfun(@times,K,dr_log'); % multiply kernel by element area
A = sparse(A); % exploit sparse structure

tools.textheader();


end



