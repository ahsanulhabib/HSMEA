function [yhat,predstd] = krig_predict(x,Surrogate)
%Function krig_predict returns the predicted response and the
%estimated prediction variance of a kriging model.

% limitations on predictor function
try
    [yhat,predvar] = krig_predictor(x,Surrogate.model);
catch
    [yhat,predvar] = krig_predictor(x,Surrogate);
end
predstd = sqrt(predvar);
return

function  [yhat,predvar] = krig_predictor(x, dmodel)
% Call:   [yhat predvar] = krig_predictor(x, dmodel)
%
% Input
% x      : trial design sites with n dimensions.
%          For mx trial sites x:
%          If mx = 1, then both a row and a column vector is accepted,
%          otherwise, x must be an mx*n matrix with the sites stored
%          rowwise.
% dmodel : Struct with DACE model; see DACEFIT
%
% Output
% yhat   : predicted response at x.
% predvar: estimated mean squared error of the kriging predictor;

% hbn@imm.dtu.dk
% Last update August 26, 2002

[m, n] = size(dmodel.S);  % number of design sites and number of dimensions
sx = size(x);            % number of trial sites and their dimension
if  min(sx) == 1 && n > 1 % Single trial point
    nx = max(sx);
    if  nx == n
        mx = 1;  x = x(:).';
    end
else
    mx = sx(1);  nx = sx(2);
end

% Normalize trial sites
x = (x - repmat(dmodel.Ssc(1,:),mx,1)) ./ repmat(dmodel.Ssc(2,:),mx,1);
q = size(dmodel.Ysc,2);  % number of response functions

if  mx == 1  % one site only
    dx = repmat(x,m,1) - dmodel.S;  % distances to design sites
    f = feval(dmodel.regr, x);
    r = feval(dmodel.corr, dmodel.theta, dx);

    % scaled dace_evaluate
    sy = f * dmodel.beta + (dmodel.gamma*r).';
    
    % prediction
    yhat = (dmodel.Ysc(1,:) + dmodel.Ysc(2,:) .* sy)';

    % prediction variance
    rt = dmodel.C \ r;
    u  = dmodel.Ft.' * rt - f.';
    v  = dmodel.G \ u;
    
    predvar = repmat(dmodel.sigma2,mx,1) .* repmat((1 + sum(v.^2) - sum(rt.^2))',1,q);

else  % several trial sites
    % get distances to design sites
    dx = zeros(mx*m,n);  kk = 1:m;
    for  k = 1 : mx
        dx(kk,:) = repmat(x(k,:),m,1) - dmodel.S;
        kk = kk + m;
    end
    
    % get regression function and correlation
    f = feval(dmodel.regr, x);
    r = feval(dmodel.corr, dmodel.theta, dx);
    r = reshape(r, m, mx);

    sy = f * dmodel.beta + (dmodel.gamma * r).';
    yhat = repmat(dmodel.Ysc(1,:),mx,1) + repmat(dmodel.Ysc(2,:),mx,1) .* sy;
    
    rt = dmodel.C \ r;
    u = dmodel.G \ (dmodel.Ft.' * rt - f.');
    predvar = repmat(dmodel.sigma2,mx,1) .* repmat((1 + colsum(u.^2) - colsum(rt.^2))',1,q);

end % of several sites

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friend functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  s = colsum(x)

if  size(x,1) == 1
    s = x;
else
    s = sum(x);
end

return
