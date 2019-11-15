
out = out_tk1;
x_norm = sqrt(sum([out.x].^2));
LAxb = [out.Axb];

for ii=1:length(out)
    Axb(ii) = norm(A*out(ii).x-b)^2;
    LAxb(ii) = norm(Lb*(A*out(ii).x-b))^2;
    Lpox(ii) = out(ii).x'*out(ii).x;
    Lprx(ii) = norm(out(ii).lambda.*out(1).Lpr*out(ii).x)^2;
    detGpo(ii) = sum(log(1./diag(out(ii).Gpo_inv)));
end

figure(11);
loglog(Axb,Lpox,'o');
text(Axb,Lpox,num2cell([out.lambda]),...
    'VerticalAlignment','bottom','HorizontalAlignment','right');

n = size(out(1).Lpr,1);
figure(12);
loglog([out.lambda],LAxb+Lprx-n.*log([out.lambda])+detGpo);
% hold on;
% loglog([out.lambda],LAxb);
% loglog([out.lambda],Lprx);
% hold off;
