function [FrontNo,MaxFNo] = NDSort(PopObj,nSort)
%NDSort - Do non-dominated sorting by efficient non-dominated sort (ENS)
%
%   FrontNo = NDSort(A,s) does non-dominated sorting on A, where A is the
%   matrix of objective values of a set of individuals, and s is the number
%   of individuals being sorted at least. FrontNo(i) means the front No. of
%   the i-th individual. The individuals have not been sorted are assigned
%   a front No. of inf.
%
%   In particular, s = 1 indicates finding only the first non-dominated
%   front, s = size(A,1)/2 indicates sorting only half of the population
%   (which is often used in the algorithm), and s = inf indicates sorting
%   the whole population.
%
%   [FrontNo,K] = NDSort(...) also returns the maximum front No. besides
%   inf.
%
%   Example:
%       [FrontNo,MaxFNo] = NDSort(PopObj,1)

%--------------------------------------------------------------------------
% The copyright of the PlatEMO belongs to the BIMK Group. You are free to
% use the PlatEMO for research purposes. All publications which use this
% platform or any code in the platform should acknowledge the use of
% "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu
% Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, 2016".
%--------------------------------------------------------------------------

% Copyright (c) 2016-2017 BIMK Group
    if nargin < 2
        nSort = size(PopObj,1);
    end
    [PopObj,~,Loc] = unique(PopObj,'rows');
    [PopObj,rank]  = sortrows(PopObj);
    Table          = hist(Loc,1:max(Loc));
    [N,M]   = size(PopObj);
    FrontNo = inf(1,N);
    MaxFNo  = 0;
    while sum(Table(FrontNo<inf)) < min(nSort,length(Loc))
        MaxFNo = MaxFNo + 1;
        for i = 1 : N
            if FrontNo(i) == inf
                Dominated = false;
                for j = i-1 : -1 : 1
                    if FrontNo(j) == MaxFNo
                        m = 2;
                        while m <= M && PopObj(i,m) >= PopObj(j,m)
                            m = m + 1;
                        end
                        Dominated = m > M;
                        if Dominated || M == 2
                            break;
                        end
                    end
                end
                if ~Dominated
                    FrontNo(i) = MaxFNo;
                end
            end
        end
    end
    FrontNo(rank) = FrontNo;
    FrontNo = FrontNo(Loc);
end