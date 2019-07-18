
% GEN_A     Generate A matrix describing kernel/transfer functions
% Author:   Timothy Sipkens, 2018-11-27
% Note:
%   Cell arrays are used for Omega_mat and Lambda_mat in order to allow for
%   the use of sparse matrices, which is necessary to store information on
%   higher resolutions grids (such as those used for phantoms).
%=========================================================================%

function A = gen_A(grid_b,grid_i,varargin)
%-------------------------------------------------------------------------%
% Inputs:
%   grid_b      Grid on which the data exists
%   grid_i      Grid on which to perform integration
%   varargin    Name-value pairs used in evaluating the CPMA tfer. func.
%-------------------------------------------------------------------------%

addpath('UBC-tfer-PMA'); % add particle mass analyzer tfer. fun. package
                         % contains dm2zp and related functions

%-- Parse measurement set points (b) -------------------------------------%
r_star = grid_b.elements;
m_star = r_star(:,1);
d_star = r_star(:,2);
n_b = grid_b.ne;
N_b = prod(n_b); % length of data vector


%-- Generate grid for intergration ---------------------------------------%
n_i = grid_i.ne;
N_i = prod(n_i); % length of integration vector

r = grid_i.elements;
m = r(:,1);
d = r(:,2);


%-- Start evaluate kernel ------------------------------------------------%
disp('Evaluating kernel...');

%-- Evaluate particle charging fractions --%
z_vec = (1:3)';
f_z = sparse(kernel.tfer_charge(d.*1e-9,z_vec)); % get fraction charged for d
n_z = length(z_vec);


%-- Evaluate DMA transfer function ---------------------------------------%
%-- Note: The DMA transfer function is 1D, speeding evaluation.
Omega_mat = cell(1,n_z); % pre-allocate for speed
for kk=1:n_z
    Omega_mat{kk} = sparse(n_b(2),n_i(2));% pre-allocate for speed
    for ii=1:n_b(2)
        Omega_mat{kk}(ii,:) = kernel.tfer_DMA(...
            grid_b.edges{2}(ii).*1e-9,...
            grid_i.edges{2}.*1e-9,...
            z_vec(kk));
    end
    
    Omega_mat{kk}(Omega_mat{kk}<(1e-7.*max(max(Omega_mat{kk})))) = 0;
        % remove numerical noise in kernel
        
	[~,jj] = max(d==grid_i.edges{2},[],2);
    Omega_mat{kk} = Omega_mat{kk}(:,jj); % extend for 2D evaluation
end


%-- Evaluate CPMA transfer function --------------------------------------%
prop_CPMA = kernel.prop_CPMA('Olfert');
disp('Evaluating CPMA contribution:');
tools.textbar(0); % initiate textbar
Lambda_mat = cell(1,n_z); % pre-allocate for speed
for kk=1:n_z
    Lambda_mat{kk} = sparse(n_b(1),N_i);% pre-allocate for speed
    for ii=1:n_b(1)
        Lambda_mat{kk}(ii,:) = kernel.tfer_CPMA(...
            grid_b.edges{1}(ii).*1e-18,m.*1e-18,...
            d.*1e-9,z_vec(kk),prop_CPMA,varargin)';
                % CPMA transfer function
        
        tools.textbar((n_b(1)*(kk-1)+ii)/(n_z*n_b(1)));
    end
end


%-- Combine to calculate kernel ------------------------------------------%
disp(' ');
disp('Compiling kernel...');
K = sparse(N_b,N_i);
for kk=1:n_z
    [~,i1] = max(m_star==grid_b.edges{1},[],2); % index corresponding to CPMA setpoint
    [~,i2] = max(d_star==grid_b.edges{2},[],2); % index correspondng to DMA setpoint
    
    K = K+f_z(z_vec(kk),:).*... % charging contribution
        Lambda_mat{kk}(i1,:).*... % CPMA contribution
        Omega_mat{kk}(i2,:); % DMA contribution
end
disp('Completed kernel evaluation.');   
disp(' ');

dr_log = grid_i.dr; % area of integral elements in [logm,logd]T space
A = bsxfun(@times,K,dr_log'); % multiply kernel by element area
A = sparse(A); % exploit sparse structure

disp('Completed A matrix calculation.');
disp(' ');


end



