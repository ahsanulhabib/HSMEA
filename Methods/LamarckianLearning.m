function [PopDecNew,PopObjNew,PopCVNew] = LamarckianLearning(Global,X0,F0,CVeps,Surro,W,Zmin,Zmax,varargin)
if Global.localimprovement
if isempty(varargin) || strcmpi(varargin{1},'S')
    SorD = 'S';
    W(W<=0) = 1e-6;
    F0_     = F0 - repmat(Zmin,size(X0,1),1);
    F0_(F0_<=0) = 1e-6;
    Zstar   = Zmin;
else
    SorD = 'D';
    W(W>=0) = -1e-6;
    F0_     = F0 - repmat(Zmax,size(X0,1),1);
    F0_(F0_>=0) = -1e-6;
    Zstar   = Zmax;
end

try
    metric = varargin{2};
catch
    metric = 'ED';
end

LB      = Global.lower;
UB      = Global.upper;
F0_uv   = F0_./repmat(sqrt(sum(F0_.^2,2)),1,size(W,2));
W_uv    = W./repmat(sqrt(sum(W.^2,2)),1,size(W,2));
Angle   = abs(acos(F0_uv*W_uv'));
[~,WId] = min(Angle,[],2);
W2W     = abs(acos(W_uv*W_uv'));
W2Wsorted = sort(W2W(WId,:),2);
UnitAngle = W2Wsorted(:,2);
if ~Global.fitcons
    [~,G0] = feval(Global.problem,Global.M,X0);
else
    G0 = pred_cons(X0,Global,Surro);
end
Options = optimoptions('fmincon','Algorithm','interior-point','display','none');
X_star = nan(size(X0));
warning('off','all');
for i = 1:size(X_star)
    V = W_uv(WId(i),:)/abs(sum(W_uv(WId(i),:),2));
    if ~isempty(G0)
        if abs(sum(max(0,G0(i,:)),2)) <= CVeps && any(~isnan(Zmax)) && all(Zmin~=Zmax)
            [X_star(i,:),~,~,output] = fmincon(@(x)obj_func(x,Global,Surro,V,Zmin,Zmax,SorD,metric),X0(i,:),[],[],[],[],LB,UB,@(x)cons_func(x,Global,Surro,V,Zstar,UnitAngle(i,:)),Options);
        else
            [X_star(i,:),~,~,output] = fmincon(@(x)violation(x,Global,Surro),X0(i,:),[],[],[],[],LB,UB,[],Options);
        end
    else
        [X_star(i,:),~,~,output] = fmincon(@(x)obj_func(x,Global,Surro,V,Zmin,Zmax,SorD,metric),X0(i,:),[],[],[],[],LB,UB,@(x)cons_func(x,Global,Surro,V,Zstar,UnitAngle(i,:)),Options);
    end
end
warning('on','all');
X_star = unique(X_star,'rows','stable');
F_star = pred_objs(X_star,Global,Surro);
Fall = [F0;F_star];
Xall = [X0;X_star];
if ~Global.fitcons
    [~,G_star] = feval(Global.problem,Global.M,X_star);
else
    G_star = pred_cons(X_star,Global,Surro);
end
PopDecNew = X_star;
PopObjNew = F_star;
PopConNew = G_star;
if ~isempty(G_star)
    PopCVNew  = abs(sum(max(0,PopConNew),2));
else
    PopCVNew  = zeros(size(F_star,1),1);
end
else
    PopDecNew = [];
    PopObjNew = [];
    PopCVNew = [];
end
end

function [obj] = obj_func(x,Global,surro,w,zmin,zmax,flag,metric)
f = pred_objs(x,Global,surro);
switch lower(metric)
    case 'aasf'
        switch flag
            case 'S'
                w(w<=0) = 1e-6;
                f_norm = (f-zmin)./(zmax-zmin);
                temp1 = f_norm./(repmat(w,size(f,1),1));
                obj = max(temp1,[],2)+(1e-3*sum(temp1,2));
            case 'D'
                w(w>=0) = -1e-6;
                f_norm = (f-zmax)./(zmax-zmin);
                temp1 = f_norm./(repmat(w,size(f,1),1));
                obj = -(max(temp1,[],2)+(1e-3*sum(temp1,2)));
        end
    case {'ed','euclidean','euclideandistance'}
        switch flag
            case 'S'
                f_norm = (f-zmin)./(zmax-zmin);
                obj = sqrt(nansum(f_norm.^2,2));
            case 'D'
                f_norm = (f-zmax)./(zmax-zmin);
                obj = -sqrt(nansum(f_norm().^2,2));
        end
end
end

function [g1,g2] = cons_func(x,Global,surro,w,z,u_angle)
f = pred_objs(x,Global,surro) - repmat(z,size(x,1),1);
if any(w < 0)
    f(f>=0) = -1e-6;
else
    f(f<=0) = 1e-6;
end
f_uv = f./repmat(sqrt(sum(f.^2,2)),1,size(w,2));
w_uv = w./repmat(sqrt(sum(w.^2,2)),1,size(w,2));
f_angle = abs(acos(f_uv*w_uv'));
c1 = (f_angle - u_angle);
if ~Global.fitcons
    [~,c2] = feval(Global.problem,Global.M,x);
else
    c2 = pred_cons(x,Global,surro);
end
g1 = [c1,c2];
g2 = [];
end

function [f] = violation(x,Global,surro)
if ~Global.fitcons
    [~,c] = feval(Global.problem,Global.M,x);
else
    c = pred_cons(x,Global,surro);
end
f = abs(sum(max(0,c),2));
end

function [y] = pred_objs(x,Global,surro)
try
    loop = numel(surro.f);
catch
    loop = numel(surro(1).f);
end
y = nan(size(x,1),loop);
for i = 1:loop
    try
        model = surro{i};
    catch
        model = surro.f{i};
    end
    modeltype = model.surro_type;
    pred_func = strcat(modeltype, '_predict');
    y(:,i) = feval(pred_func,Normalize(x,[Global.lower(:),Global.upper(:)]),model);
end
end

function [y] = pred_cons(x,Global,surro)
try
    loop = numel(surro.g);
catch
    loop = numel(surro(1).g);
end
y = nan(size(x,1),loop);
for i = 1:loop
    try
        model = surro{i};
    catch
        model = surro.g{i};
    end
    modeltype = model.surro_type;
    pred_func = strcat(modeltype, '_predict');
    y(:,i) = feval(pred_func,Normalize(x,[Global.lower(:),Global.upper(:)]),model);
end
end