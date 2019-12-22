
% TIKHONOV_BAYESF Evaluate the Bayes factor for the Tikhonov prior.
% Author: Timothy Sipkens, 2019-12-21
%=========================================================================%

function [B,F,C] = tikhonov_bayesf(A,b,x,Lpr,lambda,order)

F = -1/2*(norm(A*x-b)^2 + norm(lambda.*Lpr*x)^2);
Gpo_inv = (A')*A+lambda^2.*(Lpr')*Lpr;
C = (length(x)-order)*log(lambda) - tools.logdet(Gpo_inv)/2;
B = F+C;

end

