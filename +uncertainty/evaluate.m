
x = x_Tk1;
Gpo = Gpo_Tk1;
Gpo(Gpo<1e-3.*max(max(Gpo))) = 0;
L = L_Tk1;
Gpr = inv(L'*L);
Gb = inv(Lb'*Lb);

C = sqrt(det(2*pi.*Gpo))/(...
    sqrt(det(2*pi.*Gpr))*...
    sqrt(det(2*pi.*Gb)));

Fb = -1/2.*norm(Lb*(A*x-b));
Fpr = -1/2.*norm(L*(x-zeros(size(x))));
