function [id_fronts,f_fronts,NumComp,NumFronts] = EfficientNDSort_c(F)
if nargout == 1
    [id_fronts,~,~] = EfficientNDSort(F);
    for i = 1:numel(id_fronts)
        crowdingrank = sort_crowding(F,id_fronts{i});
        id_fronts{i} = crowdingrank(:);
    end
elseif nargout == 2
    [id_fronts,~,~] = EfficientNDSort(F);
    f_fronts = cell(numel(id_fronts),1);
    for i = 1:numel(id_fronts)
        crowdingrank = sort_crowding(F,id_fronts{i});
        id_fronts{i} = crowdingrank(:);
        f_fronts{i} = F(id_fronts{i},:);
    end
elseif nargout == 3
    [id_fronts,NumComp,~] = EfficientNDSort(F);
    f_fronts = cell(numel(id_fronts),1);
    for i = 1:numel(id_fronts)
        crowdingrank = sort_crowding(F,id_fronts{i});
        id_fronts{i} = crowdingrank(:);
        f_fronts{i} = F(id_fronts{i},:);
    end
elseif nargout == 4
    [id_fronts,NumComp,NumFronts] = EfficientNDSort(F);
    f_fronts = cell(numel(id_fronts),1);
    for i = 1:numel(id_fronts)
        crowdingrank = sort_crowding(F,id_fronts{i});
        id_fronts{i} = crowdingrank(:);
        f_fronts{i} = F(id_fronts{i},:);
    end
else
    error('ND Sort Error: Undefined Number of Output Arguments!');
end
return

function [ranks, dist] = sort_crowding(f_all, front_f)

L = length(front_f);
if L == 1
    ranks = front_f;
    dist = Inf;
else
    dist = zeros(L, 1);
    nf = size(f_all, 2);
    
    for i = 1:nf
        f = f_all(front_f, i);		% get ith objective
        [tmp, I] = sort(f);
        scale = f(I(L)) - f(I(1));
        dist(I(1)) = Inf;
        for j = 2:L-1
            id = I(j);
            id1 = front_f(I(j-1));
            id2 = front_f(I(j+1));
            if scale > 0
                dist(id) = dist(id) + (f_all(id2,i)-f_all(id1,i)) / scale;
            end
        end
    end
    dist = dist / nf;
    [tmp, I] = sort(dist, 'descend');
    ranks = front_f(I)';
end
return