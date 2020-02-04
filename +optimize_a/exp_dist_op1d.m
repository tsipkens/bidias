
% EXP_DIST_OP1D  Single parameter sensitivity study for exponential distance regularization.
%=========================================================================%

function [out] = exp_dist_op1d(A,b,d_vec,m_vec,lambda,x_ex,Gd,xi,solver)

%-- Parse inputs ---------------------------------------------%
if ~exist('solver','var'); solver = []; end
    % if computation method not specified

if ~exist('Gd','var'); Gd = []; end
if isempty(Gd); Gd = speye(2); end
     % if coordinate transform is not specified

if ~exist('xi','var'); xi = []; end % if no initial x is given
if ~exist('x_ex','var'); x_ex = []; end
%--------------------------------------------------------------%


%-- For lm/ld scaling --%
% param = logspace(log10(0.01),log10(100),42);

%-- For corr. scaling --%
param = 1-[logspace(log10(1e-3),log10(1.5),26),1.85,1.95,1.97,1.99,1.999];

disp('Optimizing exponential distance regularization:');
tools.textbar(0);
for ii=length(param):-1:1
    
    %-- For lm/ld scaling --%
    % out(ii).alpha = param(ii);
    % Gd_alt = Gd.*param(ii);
    
    %-- For corr. scaling --%
    out(ii).R12 = 1-param(ii);
    Gd_12_alt = sqrt(Gd(1,1)*Gd(2,2))*param(ii);
    Gd_alt = diag(diag(Gd))+[0,Gd_12_alt;Gd_12_alt,0];
    
    %-- Store other case parameters --%
    out(ii).lambda = lambda;
    out(ii).lm = sqrt(Gd_alt(1,1));
    out(ii).ld = sqrt(Gd_alt(2,2));
    out(ii).R12 = Gd_alt(1,2)/sqrt(Gd_alt(1,1)*Gd_alt(2,2));
    
    %-- Perform inversion --%
    [out(ii).x,~,Lpr] = invert.exp_dist(...
        A,b,d_vec,m_vec,lambda,Gd_alt,xi,solver);
    
    %-- Store ||Ax-b|| and Euclidean error --%
    if ~isempty(x_ex); out(ii).chi = norm(out(ii).x-x_ex); end
    out(ii).Axb = norm(A*out(ii).x-b);
    
    %-- Compute credence, fit, and Bayes factor --%
    [out(ii).B,out(ii).F,out(ii).C] = ...
        optimize_b.exp_dist_bayesf(A,b,out(ii).x,Lpr,lambda);
    
    tools.textbar((length(param)-ii+1)/length(param));
end


end

