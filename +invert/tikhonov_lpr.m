
% TIKHONOV_LPR Generates Tikhonov smoothing operators/matrix, L. 
% Author:   Timothy Sipkens, 2020-02-05
% 
% Inputs:
%   order       Order of the Tikhonov operator
%   n_grid      Length of first dimension of solution or Grid for x
%   x_length    Length of x vector
%               (only used if a Grid is not specified for n_grid)
%
% Outputs:
%   Lpr0        Tikhonov matrix
%=========================================================================%

function Lpr0 = tikhonov_lpr(order, n_grid, x_length)

if ~exist('order','var'); order = []; end
if isempty(order); order = 1; end

if ~isa(n_grid,'Grid')
    n = n_grid;
    if and(mod(x_length,n_grid)~=0,order~=0) % error if dimensions don't make sense
        error('Error: x_length must be integer multiple of n.');
    end
end


%-- Generate Tikhonov smoothing matrix -----------------------------------%
switch order
    case 0 % 0th order Tikhonov
        if isa(n_grid,'Grid'); Lpr0 = -speye(n_grid.Ne); % use Grid method (for partial grid support)
        else; Lpr0 = -speye(x_length);
        end
        
    case 1 % 1st order Tikhonov
        if isa(n_grid,'Grid'); Lpr0 = n_grid.l1; % use Grid method (for partial grid support)
        else
            I1 = 0.5 .* speye(n, n);
            E1 = full(sparse(1:n-1, 2:n, 1, n, n));
            D1 = E1-I1;

            m = x_length / n;
            I2 = 0.5 .* speye(m, m);
            E2 = sparse(1:m-1, 2:m, 1, m, m);
            D2 = E2 - I2;

            Lpr0 = kron(I2, D1) + kron(D2, I1);

            Lpr0 = Lpr0 - spdiags(sum(Lpr0,2), 0, x_length, x_length);
            Lpr0(end,:) = [];
        end
        
    case 1.3  % 1st order Tikhonov with rotated matrix
        slope = 3;
        if isa(n_grid,'Grid'); Lpr0 = n_grid.l1(slope); % use Grid method (for partial grid support)
        else
            I1 = 1/2 .* speye(n, n);
            E1 = full(sparse(1:n-1, 2:n, 1, n, n));
            D1 = E1-I1;
    
            m = x_length / n;
            I2 = slope/2 .* speye(m, m);
            E2 = sparse(1:m-1, 2:m, 1, m, m);
            D2 = E2 - I2;
    
            Lpr0 = kron(I2, D1) + kron(D2, I1);
    
            Lpr0 = Lpr0 - spdiags(sum(Lpr0,2), 0, x_length, x_length);
            Lpr0(end,:) = [];
        end
        
    case 2 % 2nd order Tikhonov
        if isa(n_grid,'Grid'); Lpr0 = n_grid.l2; % use Grid method (for partial grid support)
        else
            I1 = 0.25.*speye(n,n);
            E1 = sparse(1:n-1,2:n,1,n,n);
            D1 = E1+E1'-I1;

            m = x_length / n;
            I2 = 0.25.*speye(m, m);
            E2 = sparse(1:m-1, 2:m, 1, m, m);
            D2 = E2+E2'-I2;

            Lpr0 = kron(I2, D1)+kron(D2, I1);
            Lpr0 = Lpr0 - spdiags(sum(Lpr0,2), 0, x_length, x_length);
            
        end

    otherwise
        disp('The specified order of Tikhonov is not available.');
        disp(' ');
        return
end

end


