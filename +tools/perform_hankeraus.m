
% PERFORM_HANKERAUS Perform Hanke-Raus regularization parameter selection.
% Author: Timothy Sipkens, 2019-10-24
%=========================================================================%

function [lambda] = perform_hankeraus(out,A,b,bool_plot)

if ~exist('bool_plot','var'); bool_plot = []; end
if isempty(bool_plot); bool_plot = 0; end

lambda = [out.lambda];

Axb_alt = [];
for ii=length(lambda):-1:1
    Axb_alt(ii) = norm(A*out(ii).x-b);
end

if bool_plot
    figure(13);
    loglog(lambda,Axb_alt./lambda);
end

bool_continue = 0;
while ~bool_continue % loop to remove cases at end points
    [~,ind_hankeraus] = min(Axb_alt./lambda);
    bool_continue = ~(ind_hankeraus==length(lambda));
    if ~bool_continue
        lambda = lambda(1:end-1);
        Axb_alt = Axb_alt(1:end-1);
    end
end

[~,ind_hankeraus] = min(Axb_alt./lambda);
lambda = lambda(ind_hankeraus);

end

