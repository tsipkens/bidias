
% TIKHONOV_LCURVE  Automated method to find L-curve corver for Tikhonov regularization.
% Author: Arash Naseri, Timothy Sipkens, 2019-12-21\
% 
% Note: This code is written based on the methodology proposed in 
%   "A simple algorithm to find the L-curve corner" by Alessandro Cultrera 
%   and Luca Callegaro
%
% Inputs:
%   A       Model matrix
%   b       Data
%   span    Span of lambda values to consider
%   Lpr     Tikhonov matrix
%
% Ouputs:
%   x       Optimal solution
%   lambda  Optimal regularization parameter
%=========================================================================%

function [x,lambda,res_norm,x_norm] = tikhonov_lcurve(A,b,span,Lpr0)

epsilon = 1e-4; %termination threshold
phi=(1+sqrt(5))/2; % golden section

%-- Set extremes -------%
lam_log(1) = log10(min(span)); % xi=10^(lambda(i))
lam_log(4) = log10(max(span));
lam_log(2) = (lam_log(4)+phi*lam_log(1))/(phi+1);
lam_log(3) = lam_log(1)+lam_log(4)-lam_log(2);
n = (10^lam_log(4)-10^lam_log(1))/(10^lam_log(4));

for ii=1:4               
    [p(ii,1),p(ii,2),x] = ...
        l_curve_p(A,b,10^lam_log(ii),Lpr0);
end


x_norm = []; % intialize storage of solution norm
res_norm = []; % intialize storage of residual norm
while n>epsilon % main loop to refine lambda
    c2 = menger(p(1,:),p(2,:),p(3,:)); % curvature about pt. 2
    c3 = menger(p(2,:),p(3,:),p(4,:)); % curvature about pt. 3
    
    while c3<0
        lam_log(4) = lam_log(3); % contract domain towards the left
        p(4,:) = p(3,:);
        lam_log(3) = lam_log(2);
        p(3,:) = p(2,:);
        
        lam_log(2) = (lam_log(4)+phi*lam_log(1))/(phi+1);
            % determine new pt. 2
        
        [p(2,1),p(2,2),~] = ...
            l_curve_p(A,b,10^lam_log(2),Lpr0);
            % solve Tikhonov at new point
        
        c3 = menger(p(2,:),p(3,:),p(4,:));  
    end  
    
    x_norm = [x_norm;p(2,1);p(3,1);p(4,1)]; % store solution norm
    res_norm = [res_norm;p(2,2);p(3,2);p(4,2)]; % store residual norm
    
    if c2>c3
        lambda = 10^lam_log(2); % store lambda for output
        
        lam_log(4) = lam_log(3); % contract domain towards the left
        p(4,:) = p(3,:);
        lam_log(3) = lam_log(2);
        p(3,:) = p(2,:);
        
        lam_log(2) = (lam_log(4)+phi*lam_log(1))/(phi+1);
            % determine new pt. 3
        
        [p(2,1),p(2,2),x] = ...
            l_curve_p(A,b,10^lam_log(2),Lpr0);
            % only p(2,:) is recalculated
    else
        lambda = 10^lam_log(3); % store lambda for output
        
        lam_log(1) = lam_log(2); % contract domain towards the right
        p(1,:) = p(2,:);
        lam_log(2) = lam_log(3);
        p(2,:) = p(3,:);
        
        lam_log(3) = lam_log(1)+lam_log(4)-lam_log(2);
            % determine new pt. 2
        
        [p(3,1),p(3,2),x] = ...
            l_curve_p(A,b,10^lam_log(3),Lpr0);
            % only p(3,:) is recalculated
    end
    
    n = (10^lam_log(4)-10^lam_log(1))/(10^lam_log(4)); % recalculate n
end % end main loop

end



%== MENGER ===============================================================%
%   Determines the curvature of a circle by three points (Menger, 1930)
function c = menger(p1,p2,p3)

zeta(1) = p1(1,1);
eta(1) = p1(1,2);

zeta(2) = p2(1,1);
eta(2) = p2(1,2);

zeta(3) = p3(1,1);
eta(3) = p3(1,2);

pj_pk = (zeta(2)-zeta(1))^2+(eta(2)-eta(1))^2;
pk_pl = (zeta(3)-zeta(2))^2+(eta(3)-eta(2))^2;
pl_pj = (zeta(1)-zeta(3))^2+(eta(1)-eta(3))^2;

c = 2*((zeta(1)*eta(2)+zeta(2)*eta(3)+zeta(3)*eta(1)...
    -zeta(1)*eta(3)-zeta(2)*eta(1)-zeta(3)*eta(2))/...
    ((pj_pk*pk_pl*pl_pj)^.5));

end
%=========================================================================%



%== L_CURVE_P ============================================================%
function [zeta,eta,x] = l_curve_p(A,b,lambda,Lpr0)

x = invert.tikhonov(A,b,lambda,Lpr0);
zeta = norm(A*x-b,2); % residual norm
eta = norm(x(:),2); % solution norm

end
%=========================================================================%


