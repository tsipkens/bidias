
% TIKHONOV_LPR Generates Tikhonov smoothing operators/matrix, L. 
% Author:   Timothy Sipkens, 2020-04-11
% 
% Inputs:
%   order       Order of the Tikhonov operator
%   x_length    Length of x vector
%               (only used if a Grid is not specified for n_grid)
%
% Outputs:
%   Lpr0        Tikhonov matrix
%=========================================================================%

function Lpr0 = tikhonov_lpr(order,x_length)

%-- Generate Tikhonov smoothing matrix -----------------------------------%
switch order
    case 0 % 0th order Tikhonov
    	Lpr0 = speye(x_length);
    case 1 % 1st order Tikhonov
        Lpr0 = -speye(x_length);
        Lpr0 = spdiags(ones(x_length,1),1,Lpr0);
        Lpr0(end,:) = [];
    case 2 % 2nd order Tikhonov
        Lpr0 = -speye(x_length);
        Lpr0 = spdiags(0.5.*ones(x_length,2),[-1,1],Lpr0);
        Lpr0(1,:) = [];
        Lpr0(end,:) = [];
    otherwise
        disp('The specified order of Tikhonov is not available.');
        disp(' ');
        return
end

end
%=========================================================================%

