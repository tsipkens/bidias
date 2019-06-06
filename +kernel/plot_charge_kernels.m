
figure(1);

d = ([15,40,100,300].*1e-9)';
z_vec = 0:6;
T = 298;

fn_G = kernel.tfer_charge(d,z_vec,T,'Gopalakrishnan');
plot(z_vec,fn_G,'--k.');

fn_hy = kernel.tfer_charge(d,z_vec,T,'hybrid');
hold on;
plot(z_vec,fn_hy,'r');
hold off;

% hold on;
% plot(0:6,[kernel.eq_charge_frac(d(1),0,T),...
%     kernel.eq_charge_frac(d(1),1,T),...
%     kernel.eq_charge_frac(d(1),2,T),...
%     kernel.eq_charge_frac(d(1),3,T),...
%     kernel.eq_charge_frac(d(1),4,T),...
%     kernel.eq_charge_frac(d(1),5,T),...
%     kernel.eq_charge_frac(d(1),6,T)],'-o');
% hold off

fn_W = kernel.tfer_charge(d,z_vec,T,'Wiedensohler');
hold on;
plot(z_vec,fn_W,'-o');
hold off;

