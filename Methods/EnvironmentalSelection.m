function [index,SE] = EnvironmentalSelection(PopObj,PopCV,CVeps,V,Zmin,Zmax,varargin)
% The environmental selection of HSMEA
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
    if isempty(varargin) || strcmpi(varargin{1},'S')
        SorD = 'S';
        Zstar = Zmin;
    else
        SorD = 'D';
        Zstar = Zmax;
    end
    
    try
        metric = varargin{2};
    catch
        metric = 'ED';
    end
    
    [N,M] = size(PopObj);
    NV    = size(V,1);
    
    if sum(PopCV<=CVeps) > NV
        FeasIDs = find(PopCV<=CVeps);
        PopObjFeas = PopObj(FeasIDs,:);
        
        %% Translate the feasible population
        PopObjFeas = PopObjFeas - repmat(Zstar,numel(FeasIDs),1);
        if strcmpi(SorD,'S')
            PopObjFeas(PopObjFeas<=0) = 1e-6;
            PopObjFeasNorm = PopObjFeas./repmat(Zmax-Zmin,size(PopObjFeas,1),1);
        else
            PopObjFeas(PopObjFeas>=0) = -1e-6;
            PopObjFeasNorm = PopObjFeas./repmat(Zmax-Zmin,size(PopObjFeas,1),1);
        end
        
        PopObj_uv = PopObjFeas./repmat(sqrt(sum(PopObjFeas.^2,2)),[1,M]);
        V_uv      = V./repmat(sqrt(sum(V.^2,2)),[1,M]);
        
        %% Associate each solution to a reference vector
        Angle = acos(PopObj_uv*V_uv');
        [~,nonempty] = min(Angle,[],2);
        empty = setdiff((1:NV)',nonempty);
        
        %% Select one solution for each non-empty reference vector
        Next = zeros(NV,1);
        for i = unique(nonempty)'
            attached = find(nonempty==i);
            if ~isempty(attached)
                SelID    = attached(:);
                W        = V(i,:)./repmat(sqrt(sum(V(i,:).^2,2)),[1,M]);
                W        = W./repmat(abs(sum(W,2)),[1 M]);
                if strcmpi(SorD,'S')
                    W(W<=0)  = 1e-6;
                else
                    W(W>=0)  = -1e-6;
                end
                if strcmpi(metric,'AASF')
                    Temp     = PopObjFeasNorm(SelID,:)./repmat(W,numel(SelID),1);
                    Metric   = max(Temp,[],2) + 1e-3*sum(Temp,2);
                else
                    Metric   = sqrt(sum((PopObjFeasNorm(SelID,:)).^2,2));
                end
                % Select the one with the minimum AASF/Euclidean Distance value
                [~,RankMetric] = sort(Metric);
                Next(i)	= FeasIDs(SelID(RankMetric(1)));
                if strcmpi(SorD,'D')
                    Next(i) = FeasIDs(SelID(RankMetric(end)));
                end
            end
        end
        
        RemainingVecs = empty;
        while sum(Next~=0) < NV || numel(RemainingVecs) ~= 0
            [~,associated] = min(Angle(:,RemainingVecs),[],2);
            for k = unique(RemainingVecs(associated))'
                attached = find(RemainingVecs(associated)==k);
                if ~isempty(attached)
                    if strcmpi(metric,'AASF')
                        W        = V(k,:)./repmat(sqrt(sum(V(k,:).^2,2)),[1,M]);
                        W        = W./repmat(abs(sum(W,2)),[1 M]);
                        Temp     = PopObjFeasNorm(attached,:)./repmat(W,numel(attached),1);
                        Metric   = max(Temp,[],2) + 1e-3*sum(Temp,2);
                    else
                        Metric   = sqrt(sum((PopObjFeasNorm(attached,:)).^2,2));
                    end
                    %% Select the one with the best AASF/Euclidean Distance value
                    [~,RankMetric] = sort(Metric);
                    Next(k)	= FeasIDs(attached(RankMetric(1)));
                    if strcmpi(SorD,'D')
                        Next(k) = FeasIDs(attached(RankMetric(end)));
                    end
                end
            end
            RemainingVecs = setdiff(RemainingVecs,RemainingVecs(associated));
        end
    else
        [~,rank] = sort(PopCV);
        Next     = rank(1:min(N,NV));
    end
    
    % Population indices for next generation
    index = Next(Next~=0);
    SE    = SEnergy(PopObj(index,:));
end