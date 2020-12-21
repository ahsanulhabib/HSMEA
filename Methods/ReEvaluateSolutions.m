function Pop = ReEvaluateSolutions(D,Arc,Global,Model)
% Re-evaluate solutions after surrogate model update
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
Pop.decs  = D;
Pop.objs  = nan(size(D,1),Global.M);
if ~isempty(Arc.cons)
    Pop.cons = nan(size(D,1),size(Arc.cons,2));
else
    Pop.cons = [];
end
[id1,id2] = ismember(D,Arc.decs,'rows');
Pop.objs(id1,:)  = Arc.objs(id2(id2~=0),:);
if ~isempty(Pop.cons)
    Pop.cons(id1,:)  = Arc.cons(id2(id2~=0),:);
else
    Pop.cons = [];
end
if sum(~id1) > 0 && ~isempty(Model)
    for m = 1 : numel(Model.f)
        pred_func = [Model.f{m}.surro_type,'_predict'];
        Pop.objs(~id1,m) = feval(pred_func,Normalize(D(~id1,:),[Global.lower(:),Global.upper(:)]),Model.f{m});
    end
    if ~isempty(Pop.cons)
        if Global.fitcons
            for n = 1 : numel(Model.g)
                pred_func = [Model.g{n}.surro_type,'_predict'];
                Pop.cons(~id1,n) = feval(pred_func,Normalize(D(~id1,:),[Global.lower(:),Global.upper(:)]),Model.g{n});
            end
        else
            [~,C] = feval(Global.problem,Global.M,D(~id1,:));
            Pop.cons(~id1,:) = C;
        end
    end
end
end