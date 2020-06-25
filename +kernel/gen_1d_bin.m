
% GEN_1D_BIN  Evaluate equivalent kernel/transfer functions for uniformly binned data.
% Author:  Timothy Sipkens, 2020-04-21
% 
% Inputs:
%   vec_b       Vector of data bin centers
%   vec_x       Vector of reconstruction bin centers
% 
% Outputs:
%   A           Kernel matrix
%=========================================================================%

function A = gen_1d_bin(vec_b,vec_x)


%-- Start evaluate kernel ------------------------------------------------%
disp('Computing kernel...');

N_b = length(vec_b);
N_x = length(vec_x);

db = log(vec_b(2:end))-log(vec_b(1:(end-1)));
dx = log(vec_x(2:end))-log(vec_x(1:(end-1)));

nodes_b = [log(vec_b(1))-db(1)/2; ... % leading node
    (log(vec_b(2:end))+log(vec_b(1:(end-1))))./2; ... % midpoints
    log(vec_b(end))+db(end)/2]; % trailing node
nodes_x = [log(vec_x(1))-dx(1)/2; ... % leading node
    (log(vec_x(2:end))+log(vec_x(1:(end-1))))./2; ... % midpoints
    log(vec_x(end))+dx(end)/2]; % trailing node

%== Evaluate transfer function ===========================================%
K = sparse(N_b,N_x);% pre-allocate for speed
for ii=1:length(vec_b)
    K(ii,:) = max(...
        min(nodes_x(2:end),nodes_b(ii+1))-... % lower bound
        max(nodes_x(1:(end-1)),nodes_b(ii))... % upper bound
        ,0)./(nodes_b(ii+1)-nodes_b(ii)); % normalize by bin size
end

A = K.*(nodes_x(2:end)-nodes_x(1:(end-1)))'; % multiply kernel by element area (~ integration)
A = sparse(K); % exploit sparse structure

disp('Completed computing kernel matrix, <strong>A</strong>.');
disp(' ');

end



