
% GEN_KERNEL_GRID_C2  Generate A matrix describing kernel/transfer functions for PMA-SP2.
% Author: Timothy Sipkens, 2018-11-27
% 
% Notes:
% 	1. This function exploits the grid structure to minimize the number of
%	transfer function evaluations. 
%   2. Cell arrays are used for Omega_mat and Lambda_mat in order to 
%   allow for the use of sparse matrices, which is necessary to 
%   store information on higher resolutions grids
%   (such as those used for phantoms).
% 
%-------------------------------------------------------------------------%
% Inputs:
%   grid_b      Grid on which the data exists
%   grid_i      Grid on which to perform integration
%   prop_pma    Structure defining the properties of the PMA
%   varargin    Name-value pairs used in evaluating the PMA tfer. fun.
%=========================================================================%

function [A,sp] = gen_kernel_grid_c2(grid_b,grid_i,prop_pma,varargin)

if ~exist('prop_pma','var'); prop_pma = []; end
if isempty(prop_pma); prop_pma = kernel.prop_pma; end
    % import properties of PMA
    % use default properties selected by prop_pma function

    
%-- Parse measurement set points (b) -------------------------------------%
r_star = grid_b.elements;
m_star = r_star(:,2);
mrbc_star = r_star(:,1);
n_b = grid_b.ne;
N_b = prod(n_b); % length of data vector


%-- Generate grid for intergration ---------------------------------------%
n_i = grid_i.ne;
N_i = grid_i.Ne; % length of integration vector

r = grid_i.elements;
m = r(:,2);
mrbc = r(:,1);
d = (m.*1e-18./prop_pma.mass_mob_pref).^...
    (1/prop_pma.mass_mob_exp).*1e9;
    % invoke mass-mobility relation


%-- Start evaluate kernel ------------------------------------------------%
disp('Computing kernel...');

%== Evaluate particle charging fractions =================================%
z_vec = (1:3)';
f_z = sparse(kernel.tfer_charge(d.*1e-9,z_vec)); % get fraction charged for d
n_z = length(z_vec);


%== STEP 1: Evaluate DMA transfer function ===============================%
%   Note: The SP2 contribution is 1D and does not depend on the charge
%   state. It is boxcar function that takes into account discretization
%   only. 
disp('Computing SP2 contribution...');
Omega_mat = sparse(n_b(1),n_i(1));% pre-allocate for speed
for ii=1:n_b(1)
    Omega_mat(ii,:) = ...
        and(grid_b.nodes{1}(ii)<grid_i.edges{1},...
        grid_b.nodes{1}(ii+1)>grid_i.edges{1});
end

[~,jj] = max(mrbc==grid_i.edges{1},[],2);
Omega_mat = Omega_mat(:,jj);
    % repeat transfer function for repeated mass in grid_i

disp('Completed SP2 contribution.');
disp(' ');
%=========================================================================%


%== STEP 2: Evaluate PMA transfer function ===============================%
disp('Computing PMA contribution:');
tools.textbar(0); % initiate textbar
Lambda_mat = cell(1,n_z); % pre-allocate for speed
    % one cell entry per charge state
for kk=1:n_z % loop over the charge state
    Lambda_mat{kk} = sparse(n_b(1),N_i);% pre-allocate for speed
    
    for ii=1:n_b(2) % loop over m_star
        sp(ii) = tfer_pma.get_setpoint(...
            prop_pma,'m_star',grid_b.edges{2}(ii).*1e-18,varargin{:});
        Lambda_mat{kk}(ii,:) = kernel.tfer_pma(...
            sp(ii),m.*1e-18,d.*1e-9,...
            z_vec(kk),prop_pma)';
                % PMA transfer function
        
        tools.textbar((n_b(2)*(kk-1)+ii)/(n_z*n_b(2)));
    end
end
disp(' ');
%=========================================================================%


%== SETP 3: Combine to compile kernel ====================================%
disp('Compiling kernel...');
K = sparse(N_b,N_i);
for kk=1:n_z
    [~,i1] = max(m_star==grid_b.edges{2},[],2); % index corresponding to PMA setpoint
    [~,i2] = max(mrbc_star==grid_b.edges{1},[],2); % index correspondng to SP2 setpoint
    
    K = K+f_z(z_vec(kk),:).*... % charging contribution
        Lambda_mat{kk}(i1,:).*... % PMA contribution
        Omega_mat(i2,:); % SP2 contribution
end
disp('Completed kernel.');
%=========================================================================%

dr_log = grid_i.dr; % area of integral elements in [logm,logd]T space
A = bsxfun(@times,K,dr_log'); % multiply kernel by element area
A = sparse(A); % exploit sparse structure

disp('Completed computing A matrix.');
disp(' ');


end



