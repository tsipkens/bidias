
out = out_tk2;
x_norm = sqrt(sum([out.x].^2));
Axb = [out.Axb];

figure(11);
loglog(Axb,x_norm,'o');
text(Axb,x_norm,num2cell([out.lambda]),...
    'VerticalAlignment','bottom','HorizontalAlignment','right');

