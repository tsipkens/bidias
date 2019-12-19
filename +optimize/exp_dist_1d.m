
% EXP_DIST_1D  Single parameter sensitivity study for exponential distance regularization.
%=========================================================================%

function [out] = exp_dist_1d(A,b,d_vec,m_vec,lambda,x_ex,Gd,xi,solver)

%-- Parse inputs ---------------------------------------------%
if ~exist('solver','var'); solver = []; end
    % if computation method not specified

if ~exist('Gd','var'); Gd = []; end
if isempty(Gd); Gd = speye(2); end
     % if coordinate transform is not specified

if ~exist('xi','var'); xi = []; end % if no initial x is given
if ~exist('x_ex','var'); x_ex = []; end
%--------------------------------------------------------------%


% param = logspace(log10(0.01),log10(100),42); % for lm/ld scaling

param = 1-logspace(log10(1e-3),log10(1),16); % for corr. scaling
param = [param,-fliplr(param(1:(end-1)))];

disp('Optimizing exponential distance regularization:');
tools.textbar(0);
for ii=length(param):-1:1
    out(ii).param = param(ii);
    
    %-- For lm/ld scaling --%
    % Gd_alt = Gd.*param(ii);
    
    %-- For corr. scaling --%
    Gd_12_alt = sqrt(Gd(1,1)*Gd(2,2))*param(ii);
    Gd_alt = diag(diag(Gd))+[0,Gd_12_alt;Gd_12_alt,0];
    
    [out(ii).x,~,~] = invert.exp_dist(...
        A,b,d_vec,m_vec,lambda,Gd_alt,xi,solver);
    
    if ~isempty(x_ex); out(ii).chi = norm(out(ii).x-x_ex); end
    out(ii).Axb = norm(A*out(ii).x-b);
    
    tools.textbar((length(param)-ii+1)/length(param));
end


end

