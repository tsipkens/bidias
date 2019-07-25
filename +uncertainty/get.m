
% [~,~,L] = invert.tikhonov(A,b,n_x(1),1,1);
% Gb = diag(1./diag(Lb'*Lb));

% R = ones(size(L));
% s = diag(1./sqrt(diag(full(L'*L))));
% Gpr0 = s*R*s';

out = out_Tk0;
B = zeros(length(out.x(1,:)),1);
C = B; Fb = B; Fpr = B;

Gpo_fun = @(lambda) inv(lambda^2.*(full(L'*L)+eye(size(L)).*1e-2)+...
    full((Lb*A)'*Lb*A));
Gpr_fun = @(lambda) inv(lambda^2.*(full(L'*L)+eye(size(L)).*1e-2));

disp('Calcaulate the Bayes factor:');
tools.textbar(0);
for ii=1:length(B)
    
    x = out.x(:,ii);
    lambda = out.lambda(ii);
    
    Gpr = Gpr_fun(lambda);
    Gpo = Gpo_fun(lambda);
    
    eig_pr = eig(Gpr);
    eig_po = eig(Gpo);
    
    C(ii) = 1/2*sum(log(eig_po./eig_pr));
    
    Fb(ii) = -1/2.*norm(Lb*(A*x-b));
    Fpr(ii) = -1/2.*norm(L*(x-zeros(size(x))));
    
    B(ii) = Fb(ii)+Fpr(ii);
    
    tools.textbar(ii/length(B));
end

figure(50);
semilogx(out.lambda,(B-B(1))./log(10));
hold on;
semilogx(out.lambda,(Fb-Fb(1))./log(10));
semilogx(out.lambda,(Fpr-Fpr(1))./log(10));
% semilogx(out.lambda,(C-C(1))./log(10));
hold off;


