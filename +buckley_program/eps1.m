function [epsilon] = eps1(x)

epsilon = x.*erf(x) + exp(-x.^2)./(pi^0.5);