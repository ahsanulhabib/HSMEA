function [krig_model,options] = krig_train(xtr,ytr)
options = srgtsKRGSetOptions(xtr,ytr);
[krig_model,state] = srgtsKRGFit(options);
krig_model.perf = state.KRG_DACEPerf;
return

function srgtOPT = srgtsKRGSetOptions(P, T, FIT_Fn, ...
   KRG_RegressionModel, KRG_CorrelationModel, KRG_Theta0, KRG_LowerBound, KRG_UpperBound)
%Function srgtsKRGSetOptions creates the SURROGATES Toolbox option
%structure for kriging models. This structure contains the following fiels:
%
%* GENERAL PARAMETERS
%
%   SRGT   - Identifier of the surrogate technique: 'KRG'.
%   P      - NPOINTS-by-NDV matrix, where NPOINTS is the number of points
%            of the sample and NDV is the number of design variables.
%            Default: Empty matrix.
%   T      - NPOINTS-by-1 vector of responses on the P matrix points.
%            Default: Empty vector.
%   FIT_Fn - Function handle of the fitting function (which is used to
%            optimize KRG_Theta). [@dace_fit | @srgtsXVFit].
%            Default: @dace_fit.
%
%* KRIGING PARAMETERS
%
%   KRG_RegressionModel  - Function handle to a regression model. [
%                          function_handle | @dace_regpoly0 |
%                          @dace_regpoly1 | @dace_regpoly2]. Default:
%                          @dace_regpoly0.
%   KRG_CorrelationModel - Function handle to a correlation model. [
%                          function_handle | @dace_corrcubic |
%                          @dace_correxp | @dace_correxpg |
%                          @dace_corrgauss | @dace_corrlin |
%                          @dace_corrspherical | @dace_corrspline ].
%                          Default: @dace_corrgauss.
%   KRG_Theta0           - Initial guess for theta (correlation function
%                          parameters). Default:
%                          (NPOINTS^(-1/NDV))*ones(1, NDV).
%   KRG_LowerBound       - Lower bound for theta. Default: Empty vector.
%   KRG_UpperBound       - Upper bound for theta. Default: Empty vector.
%
%The SURROGATES Toolbox uses the DACE toolbox of Lophaven et al. (2002) to
%execute the kriging algorithm. As in DACE, when KRG_LowerBound and
%KRG_UpperBound are empty ([]) there will be NO optimization on theta
%(correlation function parameters).
%
%This is how you can use srgtsKRGSetOptions:
%
%     OPTIONS = srgtsKRGSetOptions: creates a structure with the empty
%     parameters.
%
%     OPTIONS = srgtsKRGSetOptions(P, T): Given the sampled data P (input
%     variables) and T (output variables), it creates a structure with
%     default parameters used for all not specified fields.
%
%     OPTIONS = srgtsKRGSetOptions(P, T, ..
%     KRG_UpperBound, FIT_Fn, FIT_LossFn KRG_RegressionModel, ...
%     KRG_CorrelationModel, KRG_Theta0, KRG_LowerBound):
%     it creates a structure with each of the specified fields.
%
%REFERENCES:
%
%Lophaven SN, Nielsen HB, and Søndergaard J, DACE - A MATLAB Kriging
%Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
%Modelling, Technical University of Denmark, 2002.
%Available at: http://www2.imm.dtu.dk/~hbn/dace/.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
srgtOPT.SRGT = 'KRG';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options
switch nargin
    case 0
        srgtOPT.P = [];
        srgtOPT.T = [];
        
        srgtOPT.FIT_Fn = [];
        
        srgtOPT.KRG_RegressionModel  = [];
        srgtOPT.KRG_CorrelationModel = [];
        srgtOPT.KRG_Theta0           = [];
        srgtOPT.KRG_LowerBound       = [];
        srgtOPT.KRG_UpperBound       = [];
        
    case 2
        [npoints,nvariables] = size(P);
        
        srgtOPT.P = P;
        srgtOPT.T = T;
        
        srgtOPT.FIT_Fn = @dace_fit;
        
        srgtOPT.KRG_RegressionModel  = @dace_regpoly0;
        srgtOPT.KRG_CorrelationModel = @dace_corrgauss;
        srgtOPT.KRG_Theta0           = (npoints^(-1/nvariables))*ones(1, nvariables);
        srgtOPT.KRG_LowerBound       = 1e-6*ones(1, nvariables);
        srgtOPT.KRG_UpperBound       = 200*ones(1, nvariables);
        
    case {3,4}
        [npoints,nvariables] = size(P);
        
        srgtOPT.P = P;
        srgtOPT.T = T;
        
        srgtOPT.FIT_Fn = @dace_fit;
        
        srgtOPT.KRG_RegressionModel  = KRG_RegressionModel;
        srgtOPT.KRG_CorrelationModel = @dace_corrgauss;
        srgtOPT.KRG_Theta0           = (npoints^(-1/nvariables))*ones(1, nvariables);
        srgtOPT.KRG_LowerBound       = 1e-6*ones(1, nvariables);
        srgtOPT.KRG_UpperBound       = 200*ones(1, nvariables);

    otherwise
        srgtOPT.P = P;
        srgtOPT.T = T;
        
        srgtOPT.FIT_Fn = FIT_Fn;
        
        srgtOPT.KRG_RegressionModel  = KRG_RegressionModel;
        srgtOPT.KRG_CorrelationModel = KRG_CorrelationModel;
        srgtOPT.KRG_Theta0           = KRG_Theta0;
        srgtOPT.KRG_LowerBound       = KRG_LowerBound;
        srgtOPT.KRG_UpperBound       = KRG_UpperBound;
end

return

function [srgtSRGT,srgtSTT] = srgtsKRGFit(srgtOPT)
%Function srgtsKRGFit fits the specified kriging model using the DACE
%toolbox of Lophaven et al. (2002).
%
%    [srgtSRGT srgtSTT] = srgtsKRGFit(srgtOPT)
%
%srgtSRGT is the surrogate structure that contains the following fields:
%* KRG_DACEModel: DACE model, a struct with the elements:
%    * regr   : function handle to the regression model.
%    * corr   : function handle to the correlation function.
%    * theta  : correlation function parameters.
%    * beta   : generalized least squares estimate.
%    * gamma  : correlation factors.
%    * sigma2 : maximum likelihood estimate of the process variance.
%    * S      : scaled design sites.
%    * Ssc    : scaling factors for design arguments.
%    * Ysc    : scaling factors for design ordinates.
%    * C      : Cholesky factor of correlation matrix.
%    * Ft     : Decorrelated regression matrix.
%    * G      : From QR factorization: Ft = Q*G'.
%
%srgtSTT is the state structure that contains the following fields:
%* FIT_Fn       : function handle of the fitting function.
%* FIT_FnVal    : value of the loss function (after fitting).
%* KRG_DACEPerf : struct with DACE performance information:
%    * nv     : Number of evaluations of objective function.
%    * perf   : (q+2)*nv array, where q is the number of elements in theta,
%               and the columns hold current values of [theta;  psi(theta);
%               type]. type = 1, 2 or 3, indicate 'start', 'explore' or
%               'move.' A negative value for type indicates an uphill step.
%
%REFERENCES:
%
%Lophaven SN, Nielsen HB, and Søndergaard J, DACE - A MATLAB Kriging
%Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
%Modelling, Technical University of Denmark, 2002.
%Available at: http://www2.imm.dtu.dk/~hbn/dace/.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch func2str(srgtOPT.FIT_Fn)
    case 'dace_fit'
        srgtSTT = srgtsFitCreateState(srgtOPT);
        if isempty(srgtOPT.KRG_LowerBound) % no optimization for theta
            [srgtSRGT.KRG_DACEModel, srgtSTT.KRG_DACEPerf, srgtSTT.FIT_FnVal] = dace_fit(...
                srgtOPT.P, ...
                srgtOPT.T, ...
                srgtOPT.KRG_RegressionModel, ...
                srgtOPT.KRG_CorrelationModel, ...
                srgtOPT.KRG_Theta0);
        else
            [srgtSRGT.KRG_DACEModel, srgtSTT.KRG_DACEPerf, srgtSTT.FIT_FnVal] = dace_fit(...
                srgtOPT.P, ...
                srgtOPT.T, ...
                srgtOPT.KRG_RegressionModel, ...
                srgtOPT.KRG_CorrelationModel, ...
                srgtOPT.KRG_Theta0, ...
                srgtOPT.KRG_LowerBound, ...
                srgtOPT.KRG_UpperBound);
        end
        
    case 'srgtsXVFit'
        [srgtSRGT, srgtSTT] = srgtsXVFit(srgtOPT);
        
end

return

function srgtSTT = srgtsFitCreateState(srgtOPT)
srgtSTT.FIT_Fn    = srgtOPT.FIT_Fn;
srgtSTT.FIT_FnVal = NaN;
return

function  [dmodel, perf, f] = dace_fit(S, Y, regr, corr, theta0, theta_lb, theta_ub)
% Constrained non-linear least-squares fit of a given correlation model to
% the provided data set and regression model
%
% Call
%   [dmodel, perf] = dace_fit(S, Y, regr, corr, theta0)
%   [dmodel, perf] = dace_fit(S, Y, regr, corr, theta0, lob, upb)
%
% Input
% S, Y    : Data points (S(i,:), Y(i,:)), i = 1,...,m
% regr    : Function handle to a regression model
% corr    : Function handle to a correlation function
% theta0  : Initial guess on theta, the correlation function parameters
% lob,upb : If present, then lower and upper bounds on theta
%           Otherwise, theta0 is used for theta
%
% Output
% dmodel  : DACE model: a struct with the elements
%    regr   : function handle to the regression model
%    corr   : function handle to the correlation function
%    theta  : correlation function parameters
%    beta   : generalized least squares estimate
%    gamma  : correlation factors
%    sigma2 : maximum likelihood estimate of the process variance
%    S      : scaled design sites
%    Ssc    : scaling factors for design arguments
%    Ysc    : scaling factors for design ordinates
%    C      : Cholesky factor of correlation matrix
%    Ft     : Decorrelated regression matrix
%    G      : From QR factorization: Ft = Q*G' .
% perf    : struct with performance information. Elements
%    nv     : Number of evaluations of objective function
%    perf   : (q+2)*nv array, where q is the number of elements
%             in theta, and the columns hold current values of
%                 [theta;  psi(theta);  type]
%             |type| = 1, 2 or 3, indicate 'start', 'explore' or 'move'
%             A negative value for type indicates an uphill step

% hbn@imm.dtu.dk
% Last update September 3, 2002

% Check design points
[m,n] = size(S);  % number of design sites and their dimension
% Check input params

if nargin == 2
    regr = @krigecr_regpoly0;
    corr = @krigecr_corrgauss;
    theta0 = (m^(-1/n))*ones(1, n);
    theta_lb = 1e-6*ones(1, n);
    theta_ub = 200*ones(1, n);
elseif nargin == 3
    corr = @krigecr_corrgauss;
    theta0 = (m^(-1/n))*ones(1, n);
    theta_lb = 1e-6*ones(1, n);
    theta_ub = 200*ones(1, n);
elseif nargin == 4
    theta0 = (m^(-1/n))*ones(1, n);
    theta_lb = 1e-6*ones(1, n);
    theta_ub = 200*ones(1, n);
elseif nargin == 5
    theta_lb = 1e-6*ones(1, n);
    theta_ub = 200*ones(1, n);
elseif nargin == 6
    theta_ub = 200*ones(1, n);
end

sY = size(Y);
if  min(sY) == 1,  
    Y = Y(:);   
    lY = max(sY);  
    sY = size(Y);
else           
    lY = sY(1); 
end
if m ~= lY
    error('S and Y must have the same number of rows');
end

% Check correlation parameters
lth = length(theta0);
if  nargin > 5  % optimization case
    if  length(theta_lb) ~= lth || length(theta_ub) ~= lth
        error('theta0, lob and upb must have the same length');
    end
    if  any(theta_lb <= 0) || any(theta_ub < theta_lb)
        error('The bounds must satisfy  0 < lob <= upb');
    end
else  % given theta
    if  any(theta0 <= 0)
        error('theta0 must be strictly positive');
    end
end

% Normalize data
mS = mean(S);   sS = std(S);
mY = mean(Y);   sY = std(Y);
% 02.08.27: Check for 'missing dimension'
j = find(sS == 0);
if  ~isempty(j),  sS(j) = 1; end
j = find(sY == 0);
if  ~isempty(j),  sY(j) = 1; end
S = (S - repmat(mS,m,1)) ./ repmat(sS,m,1);
Y = (Y - repmat(mY,m,1)) ./ repmat(sY,m,1);

% Calculate distances D between points
% mzmax = m*(m-1) / 2;        % number of non-zero distances
idvec = repmat(1:size(S,1),size(S,1),1);
idmat1 = nonzeros(tril(idvec,-1));
idmat2 = nonzeros(triu(idvec,1)');
ij = [idmat1,idmat2];
D = S(idmat1,:) - S(idmat2,:);
if  min(sum(abs(D),2) ) == 0
    SS = S; YY = Y;
% 	warning('Multiple design sites are not allowed. Selecting unique design sites...'), 
	[uS,id_S] = unique(SS,'rows','stable');
    uY = YY(id_S,:);
    [~,dist] = knnsearch(uS,uS,'k',2); 
%     nnID = id_S(dist(:,2) >= 1e-5);
    S = uS(dist(:,2) >= 1e-5,:);
    Y = uY(dist(:,2) >= 1e-5,:);
    m = size(S,1);
    % Calculate distances D between points
%     mzmax = m*(m-1) / 2;        % number of non-zero distances
    idvec = repmat(1:size(S,1),size(S,1),1);
    idmat1 = nonzeros(tril(idvec,-1));
    idmat2 = nonzeros(triu(idvec,1)');
    ij = [idmat1,idmat2];
    D = S(idmat1,:) - S(idmat2,:);
end

% Regression matrix
F = feval(regr,S);  [mF,p] = size(F);
if  mF ~= m, error('number of rows in  F  and  S  do not match'), end
if  p > mF,  error('least squares problem is underdetermined'), end

% parameters for objective function
par = struct('corr',corr, 'regr',regr, 'y',Y, 'F',F, ...
    'D', D, 'ij',ij, 'scS',sS);

% Determine theta
if  nargin > 5
    % Bound constrained non-linear optimization
    [theta, f, fit, perf] = boxmin(theta0, theta_lb, theta_ub, par);
    if  isinf(f)
        error('Bad parameter region.  Try increasing  upb'), end
else
    % Given theta
    theta = theta0(:);
    [f, fit] = objfunc(theta, par);
    perf = struct('perf',[theta; f; 1], 'nv',1);
    if  isinf(f)
        error('Bad point.  Try increasing theta0'), end
end

% Return values
dmodel = struct('regr',regr, 'corr',corr, 'theta',theta.', ...
    'beta',fit.beta, 'gamma',fit.gamma, 'sigma2',sY.^2.*fit.sigma2, ...
    'S',S, 'Ssc',[mS; sS], 'Ysc',[mY; sY], ...
    'C',fit.C, 'Ft',fit.Ft, 'G',fit.G);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friend functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [obj, fit] = objfunc(theta, par)
% Initialize
obj = inf;
fit = struct('sigma2',NaN, 'beta',NaN, 'gamma',NaN, ...
    'C',NaN, 'Ft',NaN, 'G',NaN);
m = size(par.F,1);
% Set up  R
r = feval(par.corr, theta, par.D);
idx = find(r > 0);   o = (1 : m)';
mu = (10+m)*eps;
R = sparse([par.ij(idx,1); o], [par.ij(idx,2); o], [r(idx); ones(m,1)+mu]);
% Cholesky factorization with check for pos. def.
[C, rd] = chol(R);
if  rd,  return, end % not positive definite

% Get least squares solution
C = C';   Ft = C \ par.F;
[Q, G] = qr(Ft,0);
if  rcond(G) < 1e-10
    % Check   F
    if  cond(par.F) > 1e15
        T = 'F is too ill conditioned\nPoor combination of regression model and design sites';
        error(T);
    else  % Matrix  Ft  is too ill conditioned
        return
    end
end
Yt = C \ par.y;   beta = G \ (Q'*Yt);
rho = Yt - Ft*beta;  sigma2 = sum(rho.^2)/m;
detR = prod( full(diag(C)) .^ (2/m) );
obj = sum(sigma2) * detR;
if  nargout > 1
    fit = struct('sigma2',sigma2, 'beta',beta, 'gamma',rho' / C, ...
        'C',C, 'Ft',Ft, 'G',G');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [t, f, fit, perf] = boxmin(t0, lo, up, par)
%BOXMIN  Minimize with positive box constraints

% Initialize
[t, f, fit, itpar] = start(t0, lo, up, par);
if  ~isinf(f)
    % Iterate
    p = length(t);
    if  p <= 2,  kmax = 2; else  kmax = min(p,4); end
    for  k = 1 : kmax
        th = t;
        [t, f, fit, itpar] = explore(t, f, fit, itpar, par);
        [t, f, fit, itpar] = move(th, t, f, fit, itpar, par);
    end
end
perf = struct('nv',itpar.nv, 'perf',itpar.perf(:,1:itpar.nv));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [t, f, fit, itpar] = start(t0, lo, up, par)
% Get starting point and iteration parameters

% Initialize
t = t0(:);  lo = lo(:);   up = up(:);   p = length(t);
D = 2 .^ ([1:p]'/(p+2));
ee = find(up == lo);  % Equality constraints
if  ~isempty(ee)
    D(ee) = ones(length(ee),1);   t(ee) = up(ee);
end
ng = find(t < lo | up < t);  % Free starting values
if  ~isempty(ng)
    t(ng) = (lo(ng) .* up(ng).^7).^(1/8);  % Starting point
end
ne = find(D ~= 1);

% Check starting point and initialize performance info
[f, fit] = objfunc(t,par);   nv = 1;
itpar = struct('D',D, 'ne',ne, 'lo',lo, 'up',up, ...
    'perf',zeros(p+2,200*p), 'nv',1);
itpar.perf(:,1) = [t; f; 1];
if  isinf(f)    % Bad parameter region
    return
end

if  length(ng) > 1  % Try to improve starting guess
    d0 = 16;  d1 = 2;   q = length(ng);
    th = t;   fh = f;   jdom = ng(1);
    for  k = 1 : q
        j = ng(k);    fk = fh;  tk = th;
        DD = ones(p,1);  DD(ng) = repmat(1/d1,q,1);  DD(j) = 1/d0;
        alpha = min(log(lo(ng) ./ th(ng)) ./ log(DD(ng))) / 5;
        v = DD .^ alpha;   tk = th;
        for  rept = 1 : 4
            tt = tk .* v;
            [ff, fitt] = objfunc(tt,par);  nv = nv+1;
            itpar.perf(:,nv) = [tt; ff; 1];
            if  ff <= fk
                tk = tt;  fk = ff;
                if  ff <= f
                    t = tt;  f = ff;  fit = fitt; jdom = j;
                end
            else
                itpar.perf(end,nv) = -1;   break
            end
        end
    end % improve

    % Update Delta
    if  jdom > 1
        D([1 jdom]) = D([jdom 1]);
        itpar.D = D;
    end
end % free variables

itpar.nv = nv;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [t, f, fit, itpar] = explore(t, f, fit, itpar, par)
% Explore step

nv = itpar.nv;   ne = itpar.ne;
for  k = 1 : length(ne)
    j = ne(k);   tt = t;   DD = itpar.D(j);
    if  t(j) == itpar.up(j)
        atbd = 1;   tt(j) = t(j) / sqrt(DD);
    elseif  t(j) == itpar.lo(j)
        atbd = 1;  tt(j) = t(j) * sqrt(DD);
    else
        atbd = 0;  tt(j) = min(itpar.up(j), t(j)*DD);
    end
    [ff, fitt] = objfunc(tt,par);  nv = nv+1;
    itpar.perf(:,nv) = [tt; ff; 2];
    if  ff < f
        t = tt;  f = ff;  fit = fitt;
    else
        itpar.perf(end,nv) = -2;
        if  ~atbd  % try decrease
            tt(j) = max(itpar.lo(j), t(j)/DD);
            [ff, fitt] = objfunc(tt,par);  nv = nv+1;
            itpar.perf(:,nv) = [tt; ff; 2];
            if  ff < f
                t = tt;  f = ff;  fit = fitt;
            else
                itpar.perf(end,nv) = -2;
            end
        end
    end
end % k

itpar.nv = nv;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [t, f, fit, itpar] = move(th, t, f, fit, itpar, par)
% Pattern move

nv = itpar.nv;   ne = itpar.ne;   p = length(t);
v = t ./ th;
if  all(v == 1)
    itpar.D = itpar.D([2:p 1]).^.2;
    return
end

% Proper move
rept = 1;
while  rept
    tt = min(itpar.up, max(itpar.lo, t .* v));
    [ff, fitt] = objfunc(tt,par);  nv = nv+1;
    itpar.perf(:,nv) = [tt; ff; 3];
    if  ff < f
        t = tt;  f = ff;  fit = fitt;
        v = v .^ 2;
    else
        itpar.perf(end,nv) = -3;
        rept = 0;
    end
    if  any(tt == itpar.lo | tt == itpar.up), rept = 0; end
end

itpar.nv = nv;
itpar.D = itpar.D([2:p 1]).^.25;

return
