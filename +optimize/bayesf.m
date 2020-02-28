
% BAYESF Evaluate the Bayes factor given Lpr0 and lambda.
% Author: Timothy Sipkens, 2019-12-21
%=========================================================================%

function [B,F,C] = bayesf(A,b,x,Lpr0,lambda)

r = length(x);
lx = size(A,2); % n in Thompson and Kay
lb = size(A,1); % s in Thompson and Kay

F = -1/2*(norm(A*x-b)^2 + norm(lambda.*Lpr0*x)^2); % fit

[~,~,~,S1,S2] = gsvd(full(A),full(Lpr0));
h = max(S1,[],2); % picks off diagonal elements of each row
the = diag(S2);

ct = 2*sum(log(the((lx-r+1):lb))); % contributes to the constant term
det_po = sum(log(h((lx-r+1):lb).^2 + lambda^2.*the((lx-r+1):lb).^2));

C = (r+lb-lx)*log(lambda) - det_po/2 + ct/2;
    % measurement credence
    % ignores det_b contributions

B = F+C;

end

