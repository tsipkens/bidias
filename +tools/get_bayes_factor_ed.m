
% GET_BAYES_FACTOR_ED Compute fit and approx. credence terms for Bayesian model selection.
% Applies to exponential distance prior output structure.
% Author: Timothy Sipkens, 2019-12-10
%=========================================================================%

function [Bapx,F,C] = get_bayes_factor_ed(A,b,Lb,grid_x,out)

disp('Computing fit and credence:');

tools.textbar(0);
for kk=length(out):-1:1
    
    lambda = out(kk).lambda;
    ratio = out(kk).ratio;
    R12 = out(kk).corr;
    ld = out(kk).ld;
    lm = ld*ratio;
    
    Gd = [lm^2,R12*ld^2/ratio;R12*ld^2/ratio,ld^2];
    Lpr = invert.exp_dist_lpr(grid_x.elements(:,2),...
        grid_x.elements(:,1),lambda,Gd);
    
    
    fit_b(kk) = norm(Lb*(A*out(kk).x-b))^2;
    
    fit_pr(kk) = norm(Lpr*out(kk).x)^2;
    
    cred_pr(kk) = tools.logdet(Lpr'*Lpr)/2;
    
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