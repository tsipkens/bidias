function [x,D,Lx] = exponential_distance(A,b,d_vec,m_vec,lambda,Lex,x0,solver,sigma)
% EXPONENTIAL_DISTANCE Regularization based on the exponential of the inverse distance. 
%   ...


x_length = length(A(1,:));

%-- Parse inputs ---------------------------------------------%
if ~exist('solver','var') % if computation method not specified
    solver = 'interior-point';
elseif isempty(solver)
    solver = 'interior-point';
end

if ~exist('Lex','var') % if coordinate transform is not specified
    Lex = speye(2);
elseif isempty(solver)
    Lex = speye(2);
end

if ~exist('x0','var') % if no initial x is given
    x0 = sparse(x_length,1);
elseif isempty(solver)
    x0 = sparse(x_length,1);
end
%--------------------------------------------------------------%


A = sparse(A);
b = sparse(b);

%-- Generate prior covariance matrix -----------------%
[vec_d1,vec_d2] = ndgrid(d_vec,d_vec);
[vec_m1,vec_m2] = ndgrid(m_vec,m_vec);

d1 = log(vec_m1)-log(vec_m2);
d2 = log(vec_d1)-log(vec_d2);
d = sqrt((d1.*Lex(1,1)+d2.*Lex(1,2)).^2+(d1.*Lex(2,1)+d2.*Lex(2,2)).^2); % distance

Gx = exp(-d);
if exist('sigma','var') % incorporate structure into covariance, if specified
    for ii=1:x_length
        for jj=1:x_length
            Gx(ii,jj) = Gx(ii,jj).*sigma(ii).*sigma(jj);
        end
    end
end
Gx(Gx<(0.05.*mean(mean(Gx)))) = 0; % remove any entries below thershold
Gx = Gx./max(max(Gx)); % normalize matrix structure

Gxi = inv(Gx);
Lx = chol(Gxi);
Lx = lambda.*Lx./max(max(Lx));
Lx(abs(Lx)<(0.01.*mean(mean(abs(Lx))))) = 0;
Lx = sparse(Lx);

%-- Choose and execute solver --------------------------------------------%
switch solver
    case 'interior-point' % constrained, iterative linear least squares
        options = optimoptions('lsqlin','Algorithm','interior-point','Display','none');
        x = lsqlin([A;Lx],[b;sparse(x_length,1)],...
            [],[],[],[],x0,[],[],options);
        D = []; % not specified when using this method

    case 'trust-region-reflective'
        D = ([A;Lx]'*[A;Lx])\[A;Lx]'; % invert combined matrices to get first guess
        x0 = D*[b;Lx*zeros(x_length,1)];
        
        options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
        x = lsqlin([A;Lx],[b;zeros(x_length,1)],...
            [],[],[],[],x0,[],max(x0,0),options);

    case 'algebraic' % matrix multiplication least squares (not non-negative constrained)
        D = ([A;Lx]'*[A;Lx])\[A;Lx]'; % invert combined matrices
        x = D*[b;Lx*zeros(x_length,1)];
        
    case 'algebraic-inv' % alternate algebraic least squares (less stable than previous option)
        D = inv([A;Lx]'*[A;Lx])*[A;Lx]'; % invert combined matrices
        x = D*[b;zeros(x_length,1)];
end

end

