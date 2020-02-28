
% BAYESF_B Evaluate the Bayes factor given Lpr0, lambda, and evolving A matrix.
% Author: Timothy Sipkens, 2020-02-27
%=========================================================================%

function [B,F,C] = bayesf_b(A,b,x,Lpr0,lambda)

r = length(x);
lx = size(A,2); % n in Thompson and Kay
lb = size(A,1); % s in Thompson and Kay

F = -1/2*(norm(A*x-b)^2 + norm(lambda.*Lpr0*x)^2); % fit

%-- Depreciated code ----------------------%
% Gpo_inv = (A')*A+(Lpr')*Lpr;
% C = tools.logdet(Lpr'*Lpr)/2 - tools.logdet(Gpo_inv)/2;
%------------------------------------------%

[~,~,P,S1,S2] = gsvd(full(A),full(Lpr0));
h = max(S1,[],2); % picks off diagonal elements of each row
the = diag(S2);

%-- Depreciated code ----------------------%
% ct = 2*sum(log(the((lx-r+1):lb)));
% % cp = 2*tools.logdet(P); % term cancels out between |Gpr| and |Gpo|
% det_po = sum(log(h((lx-r+1):lb).^2 + the((lx-r+1):lb).^2));
% C = det_pr/2 - det_po/2;
%------------------------------------------%

ct = 2*sum(log(the((lx-r+1):lb))); % contributes to the constant term
det_po = sum(log(h((lx-r+1):lb).^2 + lambda^2.*the((lx-r+1):lb).^2));

cp = 2*tools.logdet(P); % contributes due to presence of |Gb|
det_b = cp + 2.*sum(log(h(1:lb)));

C = (r+lb-lx)*log(lambda) - det_po/2 + ct/2 + det_b/2;
    % measurement credence, including |Gb| contribution

B = F+C;

end

