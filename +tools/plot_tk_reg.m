
out = out_tk1;
x_norm = sqrt(sum([out.x].^2));
lambda = [out.lambda];
Axb = [out.Axb];
chi0 = [out.chi];

% out = out_two_mh;
% lambda = [out.Sf];
% x_norm = sqrt(sum([out.x].^2));
% Axb = out.Axb;
% chi0 = out.chi;

Axb_alt = [];
for ii=1:length(lambda)
    Axb_alt(ii) = norm(A*out(ii).x-b);
end

%-- Plot L-curve ------------------------------------------%
figure(11);
loglog(Axb_alt,x_norm,'o');
text(Axb_alt,x_norm,num2cell(1:length(lambda)),...
    'VerticalAlignment','bottom','HorizontalAlignment','right');
hold on;
loglog(Axb_alt(2:(end-1)),...
    abs(-2.*x_norm(2:(end-1))+x_norm(3:end)+x_norm(1:(end-2))),...
    'o');
text(Axb_alt(2:(end-1)),abs(-2.*x_norm(2:(end-1))+x_norm(3:end)+x_norm(1:(end-2))),num2cell(2:(length(lambda)-1)),...
    'VerticalAlignment','bottom','HorizontalAlignment','right');
hold off;
ind_Lcurve = 20;
lambda_Lcurve = lambda(ind_Lcurve);


% figure(12);
% loglog(lambda_Tk_vec,Axb,'o');
% text(lambda_Tk_vec,Axb,num2cell(1:length(lambda_Tk_vec)),...
%     'VerticalAlignment','bottom','HorizontalAlignment','right');
% hold on;
% plot(xlim,[norm(Lb*epsilon)^2,norm(Lb*epsilon)^2]);
% hold off;
% [~,ind_Morozov] = min(abs(Axb-norm(Lb*epsilon)^2));
% lambda_Morozov = lambda_Tk_vec(ind_Morozov);

%-- Hanke-Raus rule ------------------------------------------------------%
figure(13);
loglog(lambda,Axb_alt./lambda,'o');
text(lambda,Axb_alt./lambda,num2cell(1:length(lambda)),...
    'VerticalAlignment','bottom','HorizontalAlignment','right');
ind_HankeRaus = 20; % [~,ind_HankeRaus] = min(Axb_alt./lambda);
lambda_HankeRaus = lambda(ind_HankeRaus);

figure(10);
loglog(lambda,chi0,'-');
% ylim([1e6,9e6]);
hold on;
plot([lambda_Lcurve,lambda_Lcurve],ylim);
% plot([lambda_Morozov,lambda_Morozov],ylim);
plot([lambda_HankeRaus,lambda_HankeRaus],ylim);
hold off;

%{
figure(14);
colormap(gcf,cm);
% i_vec = [11,20,27,35]; % Tk1
% i_vec = [82,94,106,118]; % Tk0
subplot(1,4,1);
bas_x.plot2d(x_st_Tk(:,i_vec(1)),0);
subplot(1,4,2);
bas_x.plot2d(x_st_Tk(:,i_vec(2)),0);
subplot(1,4,3);
bas_x.plot2d(x_st_Tk(:,i_vec(3)),0);
subplot(1,4,4);
bas_x.plot2d(x_st_Tk(:,i_vec(4)),0);
%}


