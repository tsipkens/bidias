
% GEN  Evaluate kernel/transfer functions for DMA-PMA to find A.
% Author:  Timothy Sipkens, 2020-02-04
% 
% Notes:
%   1. Cell arrays are used for Omega_mat and Lambda_mat in order to 
%   allow for the use of sparse matrices, which is necessary to 
%   store information on higher resolutions grids
%   (such as those used for phantoms).
%   2. By default, this function uses the analytical PMA transfer function 
%   corresponding to Case 1C from Sipkens et al. (Aerosol Sci. Technol. 2020b).
% 
% Inputs:
%   sp          PMA setpoint structure
%   d_star      DMA setpoints
%   grid_i      Grid on which to perform integration
%   prop_pma    Structure defining the properties of the PMA
%=========================================================================%

function A = gen(sp,d_star,grid_i,prop_pma)

if ~exist('prop_pma','var'); prop_pma = []; end
if isempty(prop_pma); prop_pma = kernel.prop_pma; end
    % import properties of PMA
    % use default properties selected by prop_pma function
if length(sp)~=length(d_star); error('Setpoint / d_star mismatch.'); end

    
%-- Parse measurement set points (b) -------------------------------------%
N_b = length(sp); % length of data vector


%-- Generate grid for intergration ---------------------------------------%
n_i = grid_i.ne;
N_i = grid_i.Ne; % length of integration vector

r = grid_i.elements;
m = r(:,1);
d = r(:,2);


%-- Start evaluate kernel ------------------------------------------------%
disp('[ Computing kernel... =============================]');

%== Evaluate particle charging fractions =================================%
z_vec = (1:3)';
f_z = sparse(kernel.tfer_charge(d.*1e-9,z_vec)); % get fraction charged for d
n_z = length(z_vec);


%== STEP 1: Evaluate DMA transfer function ===============================%
%   Note: The DMA transfer function is 1D (only a function of mobility),
%   which is exploited to speed evaluation. The results is 1 by 3 cell, 
%   with one entry per charge state.
disp('Computing DMA contribution...');
Omega_mat = cell(1,n_z); % pre-allocate for speed, one cell entry per charge state
for kk=1:n_z
    Omega_mat{kk} = sparse(N_b,n_i(2)); % pre-allocate for speed
    
    for ii=1:N_b % loop over d_star
        Omega_mat{kk}(ii,:) = kernel.tfer_dma(...
            d_star(ii).*1e-9,...
            grid_i.edges{2}.*1e-9,...
            z_vec(kk));
    end
    
    Omega_mat{kk}(Omega_mat{kk}<(1e-7.*max(max(Omega_mat{kk})))) = 0;
        % remove numerical noise in kernel
    
    [~,jj] = max(d==grid_i.edges{2},[],2);
    Omega_mat{kk} = Omega_mat{kk}(:,jj);
        % repeat transfer function for repeated mass setpoint
end
disp('Completed DMA contribution.');
disp(' ');


%== STEP 2: Evaluate PMA transfer function ===============================%
disp('Computing PMA contribution:');
tools.textbar(0); % initiate textbar
Lambda_mat = cell(1,n_z); % pre-allocate for speed
    % one cell entry per charge state
for kk=1:n_z % loop over the charge state
    Lambda_mat{kk} = sparse(N_b,N_i);% pre-allocate for speed
    
    for ii=1:N_b % loop over m_star
        Lambda_mat{kk}(ii,:) = kernel.tfer_pma(...
            sp(ii),m.*1e-18,...
            d.*1e-9,z_vec(kk),prop_pma)';
                % PMA transfer function
        
        tools.textbar((N_b*(kk-1)+ii)/(n_z*N_b));
    end
end
disp(' ');


%== SETP 3: Combine to compile kernel ====================================%
disp('Compiling kernel...');
K = sparse(N_b,N_i);
for kk=1:n_z
    K = K+f_z(z_vec(kk),:).*... % charging contribution
        Lambda_mat{kk}(:,:).*... % PMA contribution
        Omega_mat{kk}(:,:); % DMA contribution
end
disp('Kernel compiled.');

dr_log = grid_i.dr; % area of integral elements in [logm,logd]T space
A = bsxfun(@times,K,dr_log'); % multiply kernel by element area
A = sparse(A); % exploit sparse structure

disp('[ Complete ========================================]');
disp(' ');
disp(' ');


end



