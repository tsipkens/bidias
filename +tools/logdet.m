
% LOGDET  Compute the logarithm of the determinant of a matrix
function d = logdet(A)

[~,U,P] = lu(A);
du = diag(U);
c = det(P)*prod(sign(du));
d = log(c)+sum(log(abs(du)));
d = real(d);

end
