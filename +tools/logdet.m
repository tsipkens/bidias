
% LOGDET  Compute the logarithm of the determinant of a matrix
% Author: Timothy Sipkens, 2019-12-10
%=========================================================================%

function d = logdet(A)

[~,U,P] = lu(A);
du = diag(U);
c = det(P)*prod(sign(du));
d = log(c)+sum(log(abs(du)));
d = real(d);

% d = 2*sum(log(diag(chol(A))));

end
