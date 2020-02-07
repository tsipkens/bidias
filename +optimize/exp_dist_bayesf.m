
% EXP_DIST_BAYESF Evaluate the Bayes factor for the exponential distance prior.
% Author: Timothy Sipkens, 2019-12-21
%=========================================================================%

function [B,F,C] = exp_dist_bayesf(A,b,lambda,Lpr,x)

r = length(x);
lx = size(A,2); % n in Thompson and Kay
lb = size(A,1); % s in Thompson and Kay

F = -1/2*(norm(A*x-b)^2 + norm(lambda.*Lpr*x)^2);

% Gpo_inv = (A')*A+(Lpr')*Lpr;
% C = tools.logdet(Lpr'*Lpr)/2 - tools.logdet(Gpo_inv)/2;

[~,~,~,S1,S2] = gsvd(full(A),full(Lpr));
h = max(S1,[],2); % picks off diagonal elements of each row
the = diag(S2);

% ct = 2*sum(log(the((lx-r+1):lb)));
% % cp = 2*tools.logdet(P); % term cancels out between |Lpr| and |Gpo|
% det_po = sum(log(h((lx-r+1):lb).^2 + the((lx-r+1):lb).^2));
% C = det_pr/2 - det_po/2;

ct = 2*sum(log(the((lx-r+1):lb))); % contributes to the constant term
% cp = 2*tools.logdet(P); % term cancels out between |Lpr| and |Gpo|
det_po = sum(log(h((lx-r+1):lb).^2 + lambda^2.*the((lx-r+1):lb).^2));

C = (r+lb-lx)*log(lambda) - det_po/2 + ct/2;

B = F+C;

end

