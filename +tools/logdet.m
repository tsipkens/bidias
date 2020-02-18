
% LOGDET  Compute the logarithm of the determinant of a matrix
% Author: Timothy Sipkens, 2019-12-10
%=========================================================================%

function d = logdet(A)

[~,U,P] = lu(A);
u = diag(U);
c = det(P)*prod(sign(u));
d = log(c)+sum(log(abs(u)));
d = real(d);

% d = 2*sum(log(diag(chol(A))));

end
