
% EXP_DIST_BAYESF Evaluate the Bayes factor for the exponential distance prior.
% Author: Timothy Sipkens, 2019-12-21
%=========================================================================%

function [B,F,C] = exp_dist_bayesf(A,b,x,Lpr)

F = -1/2*(norm(A*x-b)^2 + norm(Lpr*x)^2);
Gpo_inv = (A')*A+(Lpr')*Lpr;
C = tools.logdet(Lpr'*Lpr)/2 - tools.logdet(Gpo_inv)/2;
B = F+C;

end

