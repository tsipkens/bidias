

out = out_Tk1;
x_norm = sum(out.x.^2);

Axb_alt = [];
for ii=1:length(out.lambda)
    Axb_alt(ii) = norm(A*out.x(:,ii)-b);
end

%-- Plot L-curve ------------------------------------------%
figure(11);
loglog(out.Axb,x_norm,'o');
text(out.Axb,x_norm,num2cell(1:length(out.lambda)),...
    'VerticalAlignment','bottom','HorizontalAlignment','right');
hold on;
loglog(out.Axb(2:(end-1)),...
    abs(-2.*x_norm(2:(end-1))+x_norm(3:end)+x_norm(1:(end-2))),...
    'o');
text(out.Axb(2:(end-1)),abs(-2.*x_norm(2:(end-1))+x_norm(3:end)+x_norm(1:(end-2))),num2cell(2:(length(out.lambda)-1)),...
    'VerticalAlignment','bottom','HorizontalAlignment','right');
hold off;
ind_Lcurve = 20;
lambda_Lcurve = out.lambda(ind_Lcurve);


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
loglog(out.lambda,Axb_alt./out.lambda,'o');
text(out.lambda,Axb_alt./out.lambda,num2cell(1:length(out.lambda)),...
    'VerticalAlignment','bottom','HorizontalAlignment','right');
[~,ind_HankeRaus] = min(Axb_alt./out.lambda);
lambda_HankeRaus = out.lambda(ind_HankeRaus);

figure(10);
loglog(out.lambda,out.chi,'.-');
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


