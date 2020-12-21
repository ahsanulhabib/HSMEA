function surrogate = TrainMultiSurrogate(Global,PopDec,PopObj,varargin)
% Train multiple surrogate models for prediction of objectives in HSMEA
%------------------------------- Reference --------------------------------
% A. Habib, H. K. Singh, T. Chugh, T. Ray and K. Miettinen, "A Multiple 
% Surrogate Assisted Decomposition-Based Evolutionary Algorithm for 
% Expensive Multi/Many-Objective Optimization," in IEEE Transactions on 
% Evolutionary Computation, vol. 23, no. 6, pp. 1000-1014, Dec. 2019, 
% DOI: 10.1109/TEVC.2019.2899030.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2016-2020 Ahsanul Habib. You are free to use HSMEA code for
% research purposes. All publications which use this code should acknowledge
% the use of "HSMEA" and reference "A. Habib, H. K. Singh, T. Chugh, T. Ray
% and K. Miettinen, "A Multiple Surrogate Assisted Decomposition-Based 
% Evolutionary Algorithm for Expensive Multi/Many-Objective Optimization,"
% in IEEE Transactions on Evolutionary Computation, vol. 23, no. 6, 
% pp. 1000-1014, Dec. 2019, DOI: 10.1109/TEVC.2019.2899030".
%--------------------------------------------------------------------------
if ~isempty(varargin)
    Theta = varargin{1};
else 
    Theta = [];
end
surrogates = Global.surrogates;
[~,D] = size(PopDec);
x_train = Normalize(PopDec,[Global.lower(:),Global.upper(:)]);
f_train = PopObj;
train_ids = (1:size(x_train,1))';
n_total = numel(train_ids);
traincount = round(n_total * Global.training_data);

if numel(train_ids) > 1.25*(2*D + 1)
    [~,~,~,~,med_ids] = kmedoids(x_train,traincount);
    t_ids = {med_ids};
    v_ids = {setdiff((1:numel(train_ids))', t_ids{:},'stable')};
else
    t_ids = {(1:numel(train_ids))'};
    v_ids = t_ids;
end

surrogate = train_model(surrogates, x_train, f_train, t_ids, v_ids, Theta);
end

function best_model = train_model(surrogates, x, f, t_ids, v_ids, Theta)
n = length(surrogates);
model_data = cell(n,1);
model_nrmse = cell(n,1);

if numel(surrogates) > 1
    for i = 1:n
        model_type = surrogates{i};
        [model_data{i}, model_nrmse{i}] = train_single_model(x, f, t_ids, v_ids, model_type, Theta);
    end
    [~,flag] = min(cell2mat(model_nrmse));
    
    best_model_nrmse = model_nrmse{flag};
    best_model_type = surrogates{flag};
    best_model = model_data{flag};
    best_model.surro_type = best_model_type;
    best_model.nrmse = best_model_nrmse;
    best_model.nrmse_all = model_nrmse;
else
    model_type = surrogates{1};
    if strcmpi(model_type,'krig')
        best_model = feval([model_type,'_train'], x, f, 'theta_0', Theta);
    else
        best_model = feval([model_type,'_train'], x, f);
    end
    best_model.surro_type = model_type;
end

end

function [model_data, model_nrmse] = train_single_model(x, f, t_ids, v_ids, model_type, Theta)
nx = size(x,2);
train_func = strcat(model_type, '_train');

if ~iscell(t_ids)
    t_ids = {t_ids};
end
if ~iscell(v_ids)
    v_ids = {v_ids};
end

id_all = []; y = []; ypred = [];
for l = 1:numel(t_ids)
    [t_ids_mod,v_ids_mod] = modify_ids(nx,t_ids{l},v_ids{l},f);
    if numel(t_ids_mod) == numel(v_ids_mod)
        y_pred = nan(numel(v_ids_mod),1);
        for i = 1:numel(v_ids_mod) % perform LOOCV
            val_id = v_ids_mod(i);
            tr_ids = setdiff(t_ids_mod,val_id,'stable');
            
            warning('off', 'all');
            if strcmp(model_type,'krig')
                model_crossval = feval(train_func, x(tr_ids,:), f(tr_ids,:), 'theta_0', Theta);
            else
                model_crossval = feval(train_func, x(tr_ids,:), f(tr_ids,:));
            end
            warning('on', 'all');
            
            y_pred(i,:) = predict_model_response(x(val_id,:), model_type, model_crossval);
        end
    else
        warning('off', 'all');
        if strcmp(model_type,'krig')
            model_crossval = feval(train_func, x(t_ids_mod,:), f(t_ids_mod,:), 'theta_0', Theta);
        else
            model_crossval = feval(train_func, x(t_ids_mod,:), f(t_ids_mod,:));
        end
        warning('on', 'all');
        
        y_pred = predict_model_response(x(v_ids_mod,:), model_type, model_crossval);
    end
    id_all = [id_all;v_ids_mod(:)];
    y = [y;f(v_ids_mod,:)];
    ypred = [ypred;y_pred];
end

if numel(id_all) > size(x,1)
    [y,ypred] = niche_y(id_all,y,ypred);
end

ydiff = y - ypred;
mse = (ydiff' * ydiff) / size(y,1);
rmse = sqrt(mse);
mm = minmax(y');
delta = mm(2) - mm(1);

if delta < 1e-6
    delta = 1;
end

model_nrmse = rmse / delta;

warning('off', 'all');
if numel(numel(t_ids)) == 1
    if numel(t_ids_mod) == numel(v_ids_mod)
        model_data = feval(train_func, x, f);
    else
        model_data = feval(train_func, x, f); % model_crossval;
    end
else
    model_data = feval(train_func, x, f);
end
warning('on', 'all');
end

function[t_ids,v_ids] = modify_ids(nx,t_ids,v_ids,y)
rejected_ids = union(find(isnan(y)),find(isinf(y)),'stable');
t_ids = setdiff(t_ids,rejected_ids,'stable');
v_ids = setdiff(v_ids,rejected_ids,'stable');
pt_required = nx + 1 - numel(t_ids);
pt_available = numel(v_ids);
if numel(t_ids) < nx+1 && pt_available > 0
    id_available = min(pt_required,pt_available);
    t_ids = [t_ids;v_ids(1:id_available,:)];
    v_ids = v_ids(id_available+1:end,:);
end
end

function [y] = predict_model_response(x, type, model)
pred_func = strcat(type, '_predict');
y = feval(pred_func, x, model);
end

function [Y,Ypred] = niche_y(id_all,y,ypred)
u_id = unique(id_all);
Y = zeros(numel(u_id),1);
Ypred = zeros(numel(u_id),1);
for i = 1:numel(u_id)
    id = id_all == u_id(i);
    Y(i) = nanmean(y(id));
    Ypred(i) = nanmean(ypred(id));
end
end