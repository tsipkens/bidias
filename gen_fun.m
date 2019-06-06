function [f,g] = gen_fun(A,b,x)

f = norm(A*x-b)^2;
g = A'*(A*x-b);

end

