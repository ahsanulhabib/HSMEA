function [Score,UsedPts] = IGD(PopObj,PF)
% <metric> <min>
% Inverted generational distance

%--------------------------------------------------------------------------
% The copyright of the PlatEMO belongs to the BIMK Group. You are free to
% use the PlatEMO for research purposes. All publications which use this
% platform or any code in the platform should acknowledge the use of
% "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu
% Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, 2016".
%--------------------------------------------------------------------------

% Copyright (c) 2016-2017 BIMK Group

    [Distance,attachment] = min(pdist2(PF,PopObj),[],2);
    Score    = mean(Distance);
    if nargout == 2
        UsedPts  = numel(unique(attachment));
    end
    if isempty(Score)
        Score = nan;
    end
end