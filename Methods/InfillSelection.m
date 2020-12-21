function [ChildSelDec,ChildRestDec,CVeps] = InfillSelection(Global,A2,Model,PopDec,PopObj,PopCV,CVeps,mu,Vi,Vn,Zmin,Zmax,varargin)
% Infill selection of HSMEA
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
    if ~isempty(varargin)
        metric = varargin{1};
    else
        metric = 'ED';
    end
    
    [ChildSelDecTi,ChildSelObjTi,ChildSelCVTi]    = SelectInfill(PopDec,PopObj,PopCV,CVeps,Vi,mu,Zmin,Zmax,'S',metric);
    [ChildSelDecTil,ChildSelObjTil,ChildSelCVTil] = LamarckianLearning(Global,ChildSelDecTi,ChildSelObjTi,CVeps,Model,Vi,Zmin,Zmax,'S',metric);
    
    if all(~isnan(Zmax)) && ~isempty(Vn)
        [ChildSelDecTn,ChildSelObjTn,ChildSelCVTn]    = SelectInfill(PopDec,PopObj,PopCV,CVeps,Vn,mu,Zmin,Zmax,'D',metric);
        [ChildSelDecTnl,ChildSelObjTnl,ChildSelCVTnl] = LamarckianLearning(Global,ChildSelDecTn,ChildSelObjTn,CVeps,Model,Vn,Zmin,Zmax,'D',metric);
    else
        ChildSelDecTn  = [];
        ChildSelObjTn  = [];
        ChildSelCVTn   = [];
        ChildSelDecTnl = [];
        ChildSelObjTnl = [];
        ChildSelCVTnl  = [];
    end
    
    ChildSelDecT        = [ChildSelDecTi;ChildSelDecTil;ChildSelDecTn;ChildSelDecTnl];
    ChildSelObjT        = [ChildSelObjTi;ChildSelObjTil;ChildSelObjTn;ChildSelObjTnl];
    ChildSelCVT         = [ChildSelCVTi;ChildSelCVTil;ChildSelCVTn;ChildSelCVTnl];
    
    [ChildSelDecT,uid]  = unique(ChildSelDecT,'rows','stable');
    ChildSelObjT        = ChildSelObjT(uid,:);
    ChildSelCVT         = ChildSelCVT(uid,:);
    
    if sum(ChildSelCVT<=CVeps) == 0
        ChildSelDec  = [];
        ChildRestDec = PopDec;
        ChildRestCV = PopCV;
        [~,rankCh]  = sort(ChildRestCV);
        [~,rankIf]  = sort(ChildSelCVT);
        numreplaced = min(mu,numel(uid));
        ChildRestDec(rankCh(numel(rankCh)-numreplaced+1:numel(rankCh)),:) = ChildSelDecT(rankIf(1:numreplaced),:);
    elseif sum(ChildSelCVT<=CVeps) <= mu && sum(ChildSelCVT<=CVeps) > 0
        ChildSelDec  = setdiff(ChildSelDecT(ChildSelCVT<=CVeps,:),A2.decs,'rows','stable');
        ChildRestDec = PopDec;
        [RestDec,restIds] = setdiff(ChildSelDecT,ChildSelDec,'rows','stable');
        RestCV       = ChildSelCVT(restIds);
        [~,rankCh]   = sort(PopCV);
        [~,rankIf]   = sort(RestCV);
        numreplaced  = min(mu,numel(restIds));
        ChildRestDec(rankCh(numel(rankCh)-numreplaced+1:numel(rankCh)),:) = RestDec(rankIf(1:numreplaced),:);
    else
        ChildDecFeas    = ChildSelDecT(ChildSelCVT<=CVeps,:);
        ChildObjFeas    = ChildSelObjT(ChildSelCVT<=CVeps,:);
        [ChildObjRank1,Rank1ID] = NDRankOne(ChildObjFeas);
        if numel(Rank1ID) <= mu
            ChildSelDec  = ChildDecFeas(Rank1ID,:);
        else
            [~,~,~,~,Mid] = kmedoids(ChildObjRank1,mu);           
            ChildSelDec  = ChildDecFeas(Rank1ID(Mid),:);
        end
        PopDec       = [PopDec;ChildSelDec];
        SelID        = knnsearch(PopDec,ChildSelDec,'k',1);
        RestID       = setdiff((1:size(PopDec,1))',SelID(:));
        ChildRestDec = PopDec(RestID,:);
    end
end