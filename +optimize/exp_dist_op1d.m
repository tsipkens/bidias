
% EXP_DIST_OP1D  Single parameter sensitivity study for exponential distance regularization.
% Optimized with respect to parameters other than lambda (which is
% is optimized by 'exp_dist_op').
% Author: Timothy Sipkens, 2019-12-20
%=========================================================================%

function [x,output] = exp_dist_op1d(A,b,lambda,Gd,grid_vec2,vec1,x_ex,xi,solver,type)

%-- Parse inputs ---------------------------------------------%
if ~exist('solver','var'); solver = []; end
    % if computation method not specified

if ~exist('Gd','var'); Gd = []; end
if isempty(Gd); Gd = speye(2); end
     % if coordinate transform is not specified

if ~exist('xi','var'); xi = []; end % if no initial x is given
if ~exist('x_ex','var'); x_ex = []; end

if ~exist('type','var'); type = []; end % if parameter to study is not given
if isempty(type); type = 'lmld'; end
%--------------------------------------------------------------%


%-- Set up parameter vector ----------%
if strcmp(type,'lmld') % lm/ld scaling
    beta_vec = logspace(log10(0.01),log10(100),42);

elseif strcmp(type,'corr') % corr. scaling
    beta_vec = 1-[logspace(log10(1e-3),log10(1.5),26),1.85,1.95,1.97,1.99,1.999];
    
else % if invalid parameter is specified
    error('Invalid parameter name.')
end
%-------------------------------------%


disp(['Optimizing exp. dist. regularization w.r.t. ',type,'...']);
tools.textbar(0);
for ii=length(beta_vec):-1:1
    
    %-- Evaluate Gd/parameters for current vector entry ------%
    if strcmp(type,'lmld') % lm/ld scaling
        output(ii).alpha = beta_vec(ii);
        Gd_alt = Gd.*beta_vec(ii);
    
    elseif strcmp(type,'corr') % corr. scaling
        output(ii).R12 = 1-beta_vec(ii);
        Gd_12_alt = sqrt(Gd(1,1)*Gd(2,2))*beta_vec(ii);
        Gd_alt = diag(diag(Gd))+[0,Gd_12_alt;Gd_12_alt,0];
    end
    %---------------------------------------------------------%
    
    %-- Store other case parameters --%
    output(ii).lambda = lambda;
    output(ii).lm = sqrt(Gd_alt(1,1));
    output(ii).ld = sqrt(Gd_alt(2,2));
    output(ii).R12 = Gd_alt(1,2)/sqrt(Gd_alt(1,1)*Gd_alt(2,2));
    
    %-- Perform inversion --%
    [output(ii).x,~,Lpr] = invert.exp_dist(...
        A,b,lambda,Gd_alt,grid_vec2,vec1,xi,solver);
    
    %-- Store ||Ax-b|| and Euclidean error --%
    if ~isempty(x_ex); output(ii).eps = norm(output(ii).x-x_ex); end
    output(ii).Axb = norm(A*output(ii).x-b);
    
    %-- Compute credence, fit, and Bayes factor --%
    [output(ii).B,output(ii).F,output(ii).C] = ...
        optimize.bayesf(A,b,output(ii).x,Lpr,lambda);
    
    tools.textbar((length(beta_vec)-ii+1)/length(beta_vec));
end


if ~isempty(x_ex)
    [~,ind_min] = min([output.eps]); % use Euclidean error
else
    [~,ind_min] = max([output.B]); % use Bayes factor
end
x = output(ind_min).x;

end

