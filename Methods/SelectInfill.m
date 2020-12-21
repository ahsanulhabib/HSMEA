function varargout = SelectInfill(PopDec,PopObj,PopCV,CVeps,V,mu,Zmin,Zmax,varargin)
% A sub-function of infill selection of HSMEA
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
        V(V<=0) = 1e-6;
    else
        SorD = 'D';
        V(V>=0) = -1e-6;
    end
    
    try 
        metric = varargin{2};
    catch
        metric = 'ED';
    end
    
    if sum(PopCV<=CVeps) > mu
        [NVa,IDva] = NoActive(PopObj(PopCV<=CVeps,:),V,Zmin,Zmax);
        
        NCluster  = min(mu,size(V,1)-NVa);
        Va        = V(IDva,:);
        
        if size(Va,1) ~= NCluster
            VaU = Va./repmat(sqrt(sum(Va.^2,2)),1,size(Va,2));
            warning off all
            IDX = kmeans(VaU,NCluster,'distance','sqeuclidean','emptyaction','singleton');
            warning on all
        else
            IDX = (1:NCluster)';
        end
        
        [N,M] = size(PopObj);
        
        if strcmpi(SorD,'S')
            PopObj = PopObj - repmat(Zmin,[N,1]);
            PopObjNorm = PopObj./repmat(Zmax-Zmin,N,1);
            PopObj(PopObj<=0) = 1e-6;
        else
            PopObj = PopObj - repmat(Zmax,[N,1]);
            PopObjNorm = PopObj./repmat(Zmax-Zmin,N,1);
            PopObj(PopObj>=0) = -1e-6;
        end
        
        PopObj_uv     = PopObj./repmat(sqrt(sum(PopObj.^2,2)),[1,M]);
        Va_uv         = Va./repmat(sqrt(sum(Va.^2,2)),[1,M]);
        Angle         = abs(acos(PopObj_uv*Va_uv'));
        [~,associate] = min(Angle,[],2);
        Cindex        = IDX(associate); % Solution to cluster
        
        Next = zeros(NCluster,1);
        
        for i = unique(Cindex)'
            current1 = find(Cindex==i & PopCV<=CVeps);
            if ~isempty(current1)
                Wu        = Va(associate(current1),:)./repmat(sqrt(sum(Va(associate(current1),:).^2,2)),[1 M]);
                [~,Wm]    = kmedoids(Wu,1);
                W         = Wm./repmat(abs(sum(Wm,2)),[1 M]);
                if strcmpi(SorD,'S')
                    W(W<=0)  = 1e-6;
                else
                    W(W>=0)  = -1e-6;
                end
                % Select the one with the minimum/maximum AASF or Euclidean Distance with respect to ideal/max point
                if strcmpi(metric,'AASF')
                    Temp      = PopObjNorm(current1,:)./repmat(W,numel(current1),1); % 
                    Metric    = max(Temp,[],2) + 1e-3*sum(Temp,2);
                else
                    Metric    = sqrt(sum((PopObjNorm(current1,:)).^2,2));
                end
                [~,best]  = min(Metric);
                if strcmpi(SorD,'D')
                    [~,best]  = max(Metric);
                end
                Next(i)   = current1(best);
            end
        end
        index     = Next(Next~=0);
    else
        [~,rank] = sort(PopCV);
        index    = rank(1:min(mu,size(PopCV,1)));
    end
    
    PopDecNew  = PopDec(index,:);
    varargout{1} = PopDecNew;
    
    if nargout >= 2
        if strcmpi(SorD,'S')
            PopObjNew = PopObj(index,:) + repmat(Zmin,[numel(index),1]);
        else
            PopObjNew = PopObj(index,:) + repmat(Zmax,[numel(index),1]);
        end
        varargout{2} = PopObjNew;
        if nargout >= 3
            PopCVNew = PopCV(index,:);
            varargout{3} = PopCVNew;
        end
    end
end