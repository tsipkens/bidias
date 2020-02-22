
% EXP_DIST_OP1D  Single parameter sensitivity study for exponential distance regularization.
%=========================================================================%

function [out] = exp_dist_op1d(A,b,lambda,Gd,d_vec,m_vec,x_ex,xi,solver,type)

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


disp('Optimizing exponential distance regularization:');
tools.textbar(0);
for ii=length(beta_vec):-1:1
    
    %-- Evaluate Gd/parameters for current vector entry ------%
    if strcmp(type,'lmld') % lm/ld scaling
        out(ii).alpha = beta_vec(ii);
        Gd_alt = Gd.*beta_vec(ii);
    
    elseif strcmp(type,'corr') % corr. scaling
        out(ii).R12 = 1-beta_vec(ii);
        Gd_12_alt = sqrt(Gd(1,1)*Gd(2,2))*beta_vec(ii);
        Gd_alt = diag(diag(Gd))+[0,Gd_12_alt;Gd_12_alt,0];
    end
    %---------------------------------------------------------%
    
    %-- Store other case parameters --%
    out(ii).lambda = lambda;
    out(ii).lm = sqrt(Gd_alt(1,1));
    out(ii).ld = sqrt(Gd_alt(2,2));
    out(ii).R12 = Gd_alt(1,2)/sqrt(Gd_alt(1,1)*Gd_alt(2,2));
    
    %-- Perform inversion --%
    [out(ii).x,~,Lpr] = invert.exp_dist(...
        A,b,lambda,Gd_alt,d_vec,m_vec,xi,solver);
    
    %-- Store ||Ax-b|| and Euclidean error --%
    if ~isempty(x_ex); out(ii).chi = norm(out(ii).x-x_ex); end
    out(ii).Axb = norm(A*out(ii).x-b);
    
    %-- Compute credence, fit, and Bayes factor --%
    [out(ii).B,out(ii).F,out(ii).C] = ...
        optimize.bayesf(A,b,out(ii).x,Lpr,lambda);
    
    tools.textbar((length(beta_vec)-ii+1)/length(beta_vec));
end


end

