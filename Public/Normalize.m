function x_norm = Normalize(x,range)
if ~isempty(x)
    N = size(x,1);
    if iscell(range)
        B = nan(size(x,2),2);
        for i = 1:numel(range)
            B(i,:) = range{i};
        end
    else
        B = range;
    end
    lb = B(:,1); LB = repmat(lb',N,1);
    ub = B(:,2); UB = repmat(ub',N,1);
    x_norm = (x-LB)./(UB-LB);
else
    x_norm = [];
end
return