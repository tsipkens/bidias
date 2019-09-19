
% TRIDIAG   Solve a tridiagonal matrix system using the Thomas algorithm.
% Author:   Timothy Sipkens, 2019-01-19
% Note:     Adapted from the Olfert laboratory.
%=========================================================================%

function x = tridiag(a,b,c,y)
%-------------------------------------------------------------------------%
% This function solves the  N x N  tridiagonal system for x:
%
%  [ b(1)  c(1)                                  ] [  x(1)  ]   [  y(1)  ] 
%  [ a(2)  b(2)  c(2)                            ] [  x(2)  ]   [  y(2)  ] 
%  [       a(3)  b(3)  c(3)                      ] [        ]   [        ] 
%  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%  [                    ...    ...    ...        ] [        ]   [        ]
%  [                        a(N-1) b(N-1) c(N-1) ] [ x(N-1) ]   [ y(N-1) ]
%  [                                a(N)   b(N)  ] [  x(N)  ]   [  y(N)  ]
%
%  y must be a vector (row or column); N is determined from its length.
%  a, b, c must be vectors of lengths N, N, and N-1 respectively.
%  a(1) is ignored.
%
%-------------------------------------------------------------------------%


%-- Check that the input arrays have acceptable sizes --------------------%
N = length(y);
if length(a)~=(N) || length(b)~=N || length(c)~=(N-1)
   error('The vectors a, b, and c must be of length N-1, N, and N-1 respectively.');
end

x = zeros(size(y)); % initialize x


%-- Phase 1: LU decomposition --------------------------------------------%
beta = b(1);
if beta==0
   error('beta = 0 at jj = 1:  matrix is singular');
end

y(1) = y(1)/b(1);

for ii = 2:N
    c(ii-1) = c(ii-1)/beta;
    beta = b(ii) - a(ii)*c(ii-1);
    y(ii) = (y(ii) - a(ii)*y(ii-1))/beta;
end


%-- Phase 2: Back-substitution -------------------------------------------%
x(N) = y(N);
for ii = (N-1):-1:1
    x(ii) = y(ii) - c(ii)*x(ii+1);
end

end
