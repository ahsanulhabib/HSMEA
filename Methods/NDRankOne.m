function [ndset,idx] = NDRankOne(set1,dir)
% Calculate Rank-1 solutions
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
% dir = 1, minimize all values
% dir = 2, maximimze all values
if nargin == 1
	dir = 1;
end
if ~isempty(set1)
    switch(dir)
        case 1
            [id,fr] = EfficientNDSort_c(set1);%dom = nd_sort_min(set1, M, N);
        case 2
            [id,fr] = EfficientNDSort_c(-set1);%dom = nd_sort_max(set1, M, N);
        otherwise
            error('wrong value of dir');
    end
    %ndset = fr{1};
    idx = id{1};

    ranks = sort_crowding(set1,idx);
    idx = ranks;
    ndset = set1(idx,:);
else
    ndset = [];
    idx = [];
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