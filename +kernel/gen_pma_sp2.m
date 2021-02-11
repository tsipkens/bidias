
% GEN_PMA_SP2  Evaluate kernel/transfer functions for PMA-SP2 to find A.
%  
%  A = kernel.gen_pma_sp2(SP,MRBC_NODES,GRID_I) evaluates the transfer 
%  function for the PMA setpoints specified by SP and the SP2 bins
%  specified by MRBC_NODES. The kernel is evaluated by integrating the 
%  transfer function over the elements in GRID_I. 
%  Please refer to the get_setpoint(...) function in tfer_pma folder for
%  more details on generating the SP struture.
% 
%  A = kernel.gen_pma_sp2(...,PROP_PMA) specifies a pre-computed PMA 
%  property data structure. If not given, the function uses the 
%  defaults of kernel.prop_pma(...).
%  
%  INPUTS:
%   sp          PMA setpoint structure
%   mrbc_nodes  Nodes (bin edges) for binned SP2 data setpoints, 
%               refractory black carbon mass
%               mrbc_edges(:,1) = lower edge;
%               mrbc_edges(:,2) = upper edge;
%   grid_i      Integration grid, e.g. grid_x
%   prop_pma    PMA properties structure, e.g. prop.r1 = 0.06
%               (Optional: default is the default properties from the
%               kernel.prop_pma function)
% 
%  OUTPUTS:
%   A           Kernel matrix
%  
%  ------------------------------------------------------------------------
%  
%  NOTES:
%  1. Cell arrays are used for Omega_mat and Lambda_mat in order to 
%  allow for the use of sparse matrices, which is necessary to 
%  store information on higher resolutions grids
%  (such as those used for phantoms).
%  2. By default, this function uses the analytical PMA transfer function 
%  corresponding to Case 1C from Sipkens et al. (Aerosol Sci. Technol. 2020b).
%  
%  AUTHOR: Timothy Sipkens, Arash Naseri, 2020-02-19

function A = gen_pma_sp2(sp, mrbc_nodes, grid_i, prop_pma)

% If not given, import default properties of PMA, 
% as selected by prop_pma function.
if ~exist('prop_pma','var'); prop_pma = []; end
if isempty(prop_pma); prop_pma = kernel.prop_pma; end
    
if length(sp)~=length(mrbc_nodes); error('Setpoint / mrbc_edges mismatch.'); end

    
%-- Parse measurement set points (b) -------------------------------------%
N_b = length(sp); % length of data vector


%-- Generate grid for intergration ---------------------------------------%
n_i = grid_i.ne;
N_i = grid_i.Ne; % length of integration vector

r = grid_i.elements;
m = r(:,2);
mrbc = r(:,1);
d = (m.*1e-18./prop_pma.rho0).^...
    (1/prop_pma.Dm).*1e9;
    % invoke mass-mobility relation


%-- Start evaluate kernel ------------------------------------------------%
tools.textheader('Computing PMA-SP2 kernel');

%== Evaluate particle charging fractions =================================%
z_vec = (1:3)';
f_z = sparse(kernel.tfer_charge(d.*1e-9,z_vec)); % get fraction charged for d
n_z = length(z_vec);



%== STEP 1: Evaluate SP2 transfer function ===============================%
%   Note: The SP2 contribution is a boxcar function that takes into 
%   account discretization only. 
disp('Computing SP2 contribution:');
tools.textbar(0); % initiate textbar
Omega_mat = sparse(N_b,n_i(1));% pre-allocate for speed
for ii=1:N_b
    Omega_mat(ii,:) = max(...
        min(grid_i.nodes{1}(2:end),mrbc_nodes(ii,2))-... % lower bound
        max(grid_i.nodes{1}(1:(end-1)),mrbc_nodes(ii,1))... % upper bound
        ,0)./(mrbc_nodes(ii,2)-mrbc_nodes(ii,1)); % normalize by SP2 bin size
    tools.textbar([ii, N_b]);
end

[~,jj] = max(mrbc==grid_i.edges{1},[],2);
Omega_mat = Omega_mat(:,jj);
    % repeat transfer function for repeated mass in grid_i

disp('Completed SP2 contribution.');
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


%== STEP 3: Combine to compile kernel ====================================%
disp('Compiling kernel...');
K = sparse(N_b,N_i);
for kk=1:n_z
    K = K+f_z(z_vec(kk),:).*... % charging contribution
        Lambda_mat{kk}(:,:).*... % PMA contribution
        Omega_mat; % SP2 contribution
end
disp('Completed kernel.');

dr_log = grid_i.dr; % area of integral elements in [logm,logd]T space
A = bsxfun(@times,K,dr_log'); % multiply kernel by element area
A = sparse(A); % exploit sparse structure

tools.textheader;

end



