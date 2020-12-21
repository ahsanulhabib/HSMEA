%% Non-dominated sorting
function [fronts,idx] = ND_Sort(f_all, id)
idx = [];
if isempty(f_all) 
	fronts = [];
	return
end

if nargin == 1
	id = (1:size(f_all,1))';
end

if isempty(id)
	fronts = [];
	return
end

try
% 	fronts = nd_sort_c(id, f_all(id,:));
    front_ids = EfficientNDSort_c(f_all(id,:));
    structure = struct('f',[]);
    fronts = repmat(structure,numel(front_ids),1);
    for f = 1:numel(front_ids)
        fronts(f).f = front_ids{f};
    end
catch
	warning('ND_SORT() MEX not available. Using slower matlab version.');
	fronts = nd_sort_m(id, f_all(id,:));
end

for i = 1:numel(fronts)
%     if i == 1
        [ranks, dist] = sort_crowding(f_all, fronts(i).f);
        idx = [idx(:);ranks(:)];
        fronts(i).f = ranks(:)';
%     else
%         idx = [idx;(fronts(i).f)'];
%     end
end

return


%% C-implementation of non-dominated sorting
function [F] = nd_sort_c(feasible, f_all)
[frontS, frontS_n] = ind_sort2(feasible', f_all');
F = [];
for i = 1:length(feasible)
	count = frontS_n(i);
	if count > 0
		tmp = frontS(1:count, i) + 1;
		F(i).f = feasible(tmp)';
	end
end
return


%% Matlab implementation of non-dominated sorting
function [F] = nd_sort_m(feasible, f_all)

front = 1;
F(front).f = [];

N = length(feasible);
M = size(f_all,2);

individual = [];
for i = 1:N
	id1 = feasible(i);
	individual(id1).N = 0;
	individual(id1).S = [];
end

% Assignging dominate flags
for i = 1:N
	id1 = feasible(i);
	f = repmat(f_all(i,:), N, 1);
	dom_less = sum(f <= f_all, 2);
	dom_more = sum(f >= f_all, 2);
	for j = 1:N
		id2 = feasible(j);
		if dom_less(j) == M && dom_more(j) < M
			individual(id1).S = [individual(id1).S id2];
		elseif dom_more(j) == M && dom_less(j) < M
			individual(id1).N = individual(id1).N + 1;
		end
	end
end

% identifying the first front
for i = 1:N
	id1 = feasible(i);
	if individual(id1).N == 0
		F(front).f = [F(front).f id1];
	end
end

% Identifying the rest of the fronts
while ~isempty(F(front).f)
	H = [];
	for i = 1 : length(F(front).f)
		p = F(front).f(i);
		if ~isempty(individual(p).S)
			for j = 1 : length(individual(p).S)
				q = individual(p).S(j);
				individual(q).N = individual(q).N - 1;
				if individual(q).N == 0
					H = [H q];
				end
			end
		end
	end
	if ~isempty(H)
		front = front+1;
		F(front).f = H;
	else
		break
	end
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