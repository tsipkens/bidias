
function A = gen_A(grid_b,grid_i)

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
f_z = kernel.tfer_charge(d.*1e-9,z_vec); % get fraction charged for d
n_z = length(z_vec);


%-- Evaluate DMA transfer function ---------------------------------------%
%-- Note: The DMA transfer function is 1D, speeding evaluation.
Omega_mat = zeros(n_b(2),n_i(2),n_z); % pre-allocate for speed
for kk=1:n_z
    for ii=1:n_b(2)
        Omega_mat(ii,:,kk) = kernel.tfer_DMA(...
            grid_b.edges{2}(ii).*1e-9,...
            grid_i.edges{2}.*1e-9,...
            z_vec(kk));
    end
end
Omega_mat(Omega_mat<(1e-8.*max(max(max(Omega_mat))))) = 0;
    % remove numerical noise in kernel
[~,jj] = max(d==grid_i.edges{2},[],2);
Omega_mat = Omega_mat(:,jj,:); % extend for 2D evaluation


%-- Evaluate CPMA transfer function --------------------------------------%
prop_CPMA = tfer_PMA.prop_CPMA('Olfert');
disp('Evaluating CPMA contribution:');
textbar(0); % initiate textbar
Lambda_mat = zeros(n_b(1),N_i,n_z); % pre-allocate for speed
for kk=1:n_z
    for ii=1:n_b(1)
        Lambda_mat(ii,:,kk) = tfer_PMA.tfer_B_diff(grid_b.edges{1}(ii).*1e-18,...
            m.*1e-18,d.*1e-9,z_vec(kk),prop_CPMA,'Rm',3)'; % CPMA transfer function
        
        textbar((n_b(1)*(kk-1)+ii)/(n_z*n_b(1)));
    end
end


%-- Combine to calculate kernel ------------------------------------------%
disp('Compiling kernel...');
K = zeros(N_b,N_i);
for kk=1:n_z
    [~,i1] = max(m_star==grid_b.edges{1},[],2);
    [~,i2] = max(d_star==grid_b.edges{2},[],2);
    K = K+f_z(z_vec(kk),:).*... % charging contribution
        Lambda_mat(i1,:,kk).*... % CPMA contribution
        Omega_mat(i2,:,kk); % DMA contribution
end
disp('Completed kernel evaluation.');   
disp(' ');

dr_log = grid_i.dr; % area of integral elements in [lnm,lnd] space
A = bsxfun(@times,K,dr_log'); % multiply kernel by element area
A = sparse(A); % exploit sparse structure

disp('Completed A matrix calculation.');
disp(' ');


end



