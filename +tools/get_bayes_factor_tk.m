
% GET_BAYES_FACTOR_TK Compute fit and approx. credence terms for Bayesian model selection.
% Also outputs index of optimal reconstruction based on this scheme.
% Applies to Tikhonov regularization output structure.
% Author: Timothy Sipkens, 2019-12-10
%=========================================================================%

function [Bapx,F,C] = get_bayes_factor_tk(A,b,Lb,out)

disp('Computing fit and credence:');

cred_pr = length(out(1).x).*log([out.lambda]);
tools.textbar(0);
for kk=length(out):-1:1
    fit_b(kk) = norm(Lb*(A*out(kk).x-b))^2;
    
    Lpr = out(kk).lambda.*out(1).Lpr;
    fit_pr(kk) = norm(Lpr*out(kk).x)^2;
    
    Gpo_inv = (Lb*A)'*(Lb*A)+(Lpr')*Lpr;
    cred_po(kk) = tools.logdet(Gpo_inv)/2;
    
    tools.textbar((length(out)-kk+1)./length(out));
end

F = -1/2.*(fit_b+fit_pr);
C = cred_pr-cred_po;
Bapx = F+C;

disp('Complete.');
disp(' ');

end