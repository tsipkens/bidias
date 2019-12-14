
% GET_POSTERIOR Evaluates the posterior covariance matrix. 
% This includes the posterior std. dev. for each point in the reconstruction.
% Lpr should include lambda.
% Author: Timothy Sipkens, 2019-12-09
%=========================================================================%

function [Gpo,spo] = get_posterior(A,Lb,Lpr)

Gpo_inv = (Lb*A)'*(Lb*A)+(Lpr')*Lpr;
Gpo = inv(Gpo_inv);
spo = sqrt(diag(Gpo));

end


