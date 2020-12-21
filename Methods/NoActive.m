function [num,active] = NoActive(PopObj,V,varargin)
% Detect inactive reference vectors
%------------------------------- Copyright --------------------------------
% Copyright (c) 2016-2017 BIMK Group.
%------------------------------- Reference --------------------------------
% "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    if ~isempty(varargin{1})
        Zmin = varargin{1};
    else
        Zmin = min(PopObj,[],1);
    end
    
    if ~isempty(varargin{2})
        Zmax = varargin{2};
    else
        Zmax = max(PopObj,[],1);
    end
        
    [N,M] = size(PopObj);
    NV    = size(V,1);

    %% Translate the population
    if sum(V(:,1)<0) == NV
        PopObj = PopObj - repmat(Zmax,N,1);
        PopObj(PopObj>=0) = -1e-6;
    else
        PopObj = PopObj - repmat(Zmin,N,1);
        PopObj(PopObj<=0) = 1e-6;
    end
    %% Calculate corresponding unit vectors
    PopObj_uv = PopObj./repmat(sqrt(sum(PopObj.^2,2)),[1,M]);
    V_uv      = V./repmat(sqrt(sum(V.^2,2)),[1,M]);
    %% Associate each solution to a reference vector
    %Angle   = acos(1-pdist2(PopObj,V,'cosine'));
    Angle   = acos(PopObj_uv*V_uv');
    [~,associate] = min(Angle,[],2);
    active  = unique(associate);
	num     = NV-length(active);
end