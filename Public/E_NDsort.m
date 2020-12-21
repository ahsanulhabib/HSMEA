function [id_fronts_all,F_fronts_all] = E_NDsort(F)
if ~isempty(F)
    fronts = EfficientNDSort_c(F);
    for i = 1:numel(fronts)
        id_fronts_all{i,1} = fronts{i};
        F_fronts_all{i,1} = F(fronts{i},:);
    end
else
    id_fronts_all{1} = [];
    F_fronts_all{1} = [];
end
end