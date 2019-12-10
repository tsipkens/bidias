
% GET_POSTERIOR Evaluates the posterior covariance matrix 
% This includes the std. dev. for each point.
% Author: Timothy Sipkens, 2019-09-12
%=========================================================================%

function [Gpo,spo] = get_posterior(A,Lb,Lpr)

disp('Computing posterior covariance matrix...');

Gpo_inv = (Lb*A)'*(Lb*A)+(Lpr')*Lpr;
Gpo = inv(Gpo_inv);
spo = sqrt(diag(Gpo));

disp('Complete.');
disp(' ');

end


