
% MARKOWSKI Funtion to perform Markowski/Buckley-type smoothing.
%   Used in the Twomey-Markowski algorithm.
% Author: Timothy Sipkens, 2020-02-06
% Note: n can be an integer or a Grid
%=========================================================================%

function [x,G_smooth,SIGMA] = markowski(A,b,Lb,x,n,iter,opt_smooth,Sf,SIGMA_end)

x_length = length(x);

if strcmp('Buckley',opt_smooth) % smoothing according to Buckley
    G_smooth = G_Buckley(n,x_length,Sf);
elseif strcmp('Grid',opt_smooth) % smoothing based on adjacency matrix in Grid
    G_smooth = G_grid(n,x_length,Sf);
else % original Markowski-type smoothing matrix
    G_smooth = G_Markowski(n,x_length);
end

%-- Perform smoothing over multiple iterations ---------------------------%
jj = 1;
while jj<=iter
    x = G_smooth*x; % apply smoothing

    SIGMA = calc_mean_sq_error(Lb*A,x,Lb*b); % calculate mean square error
    if SIGMA>SIGMA_end % exit smoothing if mean square error exceeds end position
        disp(['SMOOTHING: Completed smoothing algorithm after ',num2str(jj),' iteration(s).']);
        return
    end

    jj = jj+1;
end

disp(['SMOOTHING algorithm did not necessarily converge after ',num2str(iter),' iteration(s).']);

end



%== G_MARKOWSKI ==========================================================%
%   Generates a smoothing matrix of the form based on original work of Markowski
function G_smooth = G_Markowski(n,x_length)

G_smooth = 0.5.*eye(x_length);
for jj=1:x_length
    if ~(mod(jj,n)==0)
        G_smooth(jj,jj+1) = 0.125;
    else
        G_smooth(jj,jj) = G_smooth(jj,jj)+0.125;
    end

    if ~(mod(jj-1,n)==0)
        G_smooth(jj,jj-1) = 0.125;
    else
        G_smooth(jj,jj) = G_smooth(jj,jj)+0.125;
    end

    if jj<=(x_length-n)
        G_smooth(jj,jj+n) = 0.125;
    else
        G_smooth(jj,jj) = G_smooth(jj,jj)+0.125;
    end

    if jj>n
        G_smooth(jj,jj-n) = 0.125;
    else
        G_smooth(jj,jj) = G_smooth(jj,jj)+0.125;
    end
end

end
%=========================================================================%


%== G_BUCKLEY ============================================================%
%   Generates the smoothing matrix suggested by Buckley et al. (2017)
function G_smooth = G_Buckley(n,x_length,Sf)

norm_factor = 0.5+8*Sf;
G_smooth = 0.5/norm_factor.*eye(x_length);
for jj=1:x_length
    if mod(jj,n)==0 % right side
        if jj>(x_length-n) % bottom, right corner
            G_smooth(jj,jj) = 0.75;
            G_smooth(jj,jj-1) = 0.125;
            G_smooth(jj,jj-n) = 0.125;
        elseif jj<=n % top, right corner
            G_smooth(jj,jj) = 0.75;
            G_smooth(jj,jj-1) = 0.125;
            G_smooth(jj,jj+n) = 0.125;
        else % in the middle, right
            G_smooth(jj,jj-1) = 0.25;
            G_smooth(jj,jj) = 0.75;
        end
    elseif mod(jj-1,n)==0  % left side
        if jj>(x_length-n) % bottom, left corner
            G_smooth(jj,jj) = 0.75;
            G_smooth(jj,jj+1) = 0.125;
            G_smooth(jj,jj-n) = 0.125;
        elseif jj<=n % top, left corner
            G_smooth(jj,jj) = 0.75;
            G_smooth(jj,jj+1) = 0.125;
            G_smooth(jj,jj+n) = 0.125;
        else % in the middle, left
            G_smooth(jj,jj+1) = 0.25;
            G_smooth(jj,jj) = 0.75;
        end
    else % middle
        if jj>(x_length-n) % bottom, middle
            G_smooth(jj,jj-n) = 0.25;
            G_smooth(jj,jj) = 0.75;
        elseif jj<=n % top, middle
            G_smooth(jj,jj+n) = 0.25;
            G_smooth(jj,jj) = 0.75;
        else % in the middle, middle
            G_smooth(jj,jj-1) = Sf/norm_factor;
            G_smooth(jj,jj+1) = Sf/norm_factor;
            G_smooth(jj,jj+n) = Sf/norm_factor;
            G_smooth(jj,jj-n) = Sf/norm_factor;
            G_smooth(jj,jj+1+n) = Sf/norm_factor;
            G_smooth(jj,jj+1-n) = Sf/norm_factor;
            G_smooth(jj,jj-1+n) = Sf/norm_factor;
            G_smooth(jj,jj-1-n) = Sf/norm_factor;
        end
    end
end

end
%=========================================================================%


%== G_GRID ===============================================================%
%   Generates the smoothing matrix when a partial grid is provided.
%   Uses a four-point stencil. 
function G_smooth = G_grid(grid,~,Sf)

g1 = grid.adj; % off-diagonal components
norm_factor = 0.5+sum(g1,2).*Sf;

G_smooth = (0.5.*speye(size(g1))+Sf.*g1)./norm_factor;

end
%=========================================================================%



%== CALC_MEAN_SQ_ERR =====================================================%
%   Calculates the mean squared error of the non-zero entries of b
function SIGMA = calc_mean_sq_error(A,x,b)

sqErr = (A*x-b).^2; % squared errors, normalized by expected error in b
SIGMA = mean(sqErr(b~=0)); % average square error for cases where b~= 0

end
%=========================================================================%


