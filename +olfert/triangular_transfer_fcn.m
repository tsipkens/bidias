function [Omega]=triangular_transfer_fcn(beta,x_tilde)


Omega = 1/(2*beta)*(abs(x_tilde-(1+beta)) + ...
    abs(x_tilde-(1-beta)) - 2*abs(x_tilde-1));




