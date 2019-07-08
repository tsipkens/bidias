
function [x_t,grid_t,mg] = gen_phantom(param,span)
% GEN_PHANTOM Generates a mass-mobiltiy distribution phantom.
% Author:   Timothy Sipkens, 2018-12-04

n_t = [540,550]; % resolution of phantom distribution
grid_t = Grid(span,... 
    n_t,'logarithmic'); % generate grid of which to represent phantom

[x_t,mg] = p_fun(grid_t,param);

end


function [p,mg] = p_fun(grid_t,param)

%-- Parse inputs ---------------------------------------------------------%
for ll=1:length(param) % loop over different modes
    if ~isfield(param(ll),'opt_m')
        param(ll).opt_m = 'logn';
    elseif isempty(param(ll).opt_m)
        param(ll).opt_m = 'logn';
    end
end


%-- Evaluate phantom mass-mobility distribution --------------------------%
m0 = grid_t.elements(:,1);
d0 = grid_t.elements(:,2);

rho = @(d,k,Dm) 6*k./(pi*d.^(3-Dm));

mg = @(d0,ll) log(1e-9.*rho(d0,param(ll).k,param(ll).Dm).*...
    pi/6.*(d0.^3)); % geometric mean mass in fg

p = zeros(length(m0),1);
for ll=1:length(param) % loop over different modes
    if strcmp(param(ll).opt_m,'norm')
        p_m = normpdf(m0,...
            exp(mg(d0,ll)),param(ll).sm.*exp(mg(d0,ll)));
    else
        p_m = lognpdf(m0,mg(d0,ll),log(param(ll).sm));
    end
    
    p_temp = lognpdf(d0,log(param(ll).dg),log(param(ll).sg)).*...
        p_m;
    p = p+p_temp;
end

p = p./length(param);
p = p.*(d0.*m0); % convert to [lnm,lnd]T space

end


