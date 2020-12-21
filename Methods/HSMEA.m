function HSMEA(Global)
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
    rng(Global.baseseed+Global.run,'twister');
    mu = Global.mu;
    if ~Global.localimprovement
        Global.T = 50;
    else
        Global.T = 1;
    end
    T = Global.T;
    tic;
    
    %% Generate the reference points and population
    if Global.M == 2
        V0 = DirectionVector(Global.M,[29 0]);
        Global.N = size(V0,1);
    elseif Global.M == 3
        V0 = DirectionVector(Global.M,[13 0]);
        Global.N = size(V0,1);
    elseif Global.M == 4
        V0 = DirectionVector(Global.M,[7 0]);
        Global.N = size(V0,1);
    elseif Global.M == 5
        V0 = DirectionVector(Global.M,[5 0]);
        Global.N = size(V0,1);
    elseif Global.M == 6
        V0 = DirectionVector(Global.M,[4 1]);
        Global.N = size(V0,1);
    elseif Global.M >= 7
        V0 = DirectionVector(Global.M,[3 2]);
        Global.N = size(V0,1);
    else
        [V0,Global.N] = UniformPoint(Global.N,Global.M);
    end
    
	Vi       = V0;
    NI       = 11*Global.D-1;
    P        = repmat(Global.lower,NI,1) + repmat((Global.upper-Global.lower),NI,1).*lhsdesign(NI,Global.D);
    [~,C]    = feval(Global.problem,Global.M,P(1,:));
    
    if ~isempty(C)
        NI    = 50; clear A2;
        A2    = EvaluateSolutions(P,Global);
    else
        A2    = EvaluateSolutions(P,Global);
    end
    
    A1.decs = A2.decs;
    A1.objs = A2.objs;
    A1.cons = A2.cons;
    
    A1Obj = A1.objs;
    A1Con = A1.cons;
    if isempty(A1Con)
        A1CV  = zeros(size(A1Obj,1),1);
    else
        A1CV  = abs(sum(max(0,A1Con),2));
    end
    
    CVeps = mean(A1CV)*(sum(A1CV==0)/numel(A1CV));
    
    Zmax  = nan(1,Global.M);
    if sum(A1CV<=CVeps) > 1
        Zmin  = min(A1Obj(A1CV<=CVeps,:),[],1);
        Zmax  = max(A1Obj(A1CV<=CVeps,:),[],1);
        if all(Zmin ~= Zmax)
            Vi(1:Global.N,:) = V0.*repmat(Zmax-Zmin,size(V0,1),1);
        end
    else
        Zmin  = min(A1Obj,[],1);
        Vi(1:Global.N,:) = V0;
    end
    if Global.numRVset > 1
        Vn = -Vi;
    else
        Vn = [];
    end

    THETAobj = 5*ones(Global.M,Global.D);
    Model.f = cell(1,Global.M);

    f_train	= A2.objs;
    for m = 1 : Global.M
        dmodel     = TrainMultiSurrogate(Global,A2.decs,f_train(:,m),THETAobj(m,:));
        Model.f{m} = dmodel;
        if strcmpi(dmodel.surro_type,'krig')
            THETAobj(m,:) = dmodel.theta;
        end
    end
    
    if Global.fitcons
        THETAcon = 5*ones(size(A2.cons,2),Global.D);
        Model.g = cell(1,size(A2.cons,2));
        g_train	= A2.cons;
        for n = 1 : size(A2.cons,2)
            dmodel     = TrainMultiSurrogate(Global,A2.decs,g_train(:,n),THETAcon(n,:));
            Model.g{n} = dmodel;
            if strcmpi(dmodel.surro_type,'krig')
                THETAcon(n,:) = dmodel.theta;
            end
        end
    end
    
    Gen = 1; NumOffsprings = [];
    %% Optimization Loop
    while size(A2.decs,1) < Global.FEmax
        done = round(((size(A2.decs,1)/Global.FEmax)*100)*100)/100;
        disp([Global.problem,':  ',num2str(Global.M),'obj. ','Run-',num2str(Global.run),' Gen-',num2str(Gen),' (',num2str(done),' % finished)' ]);
        % Refresh the model and generate promising solutions
        A1Dec  = A1.decs;
        A1Obj  = A1.objs;
        A1Con  = A1.cons;
        PopDec = A1Dec;
        PopObj = A1Obj;
        PopCon = A1Con;
        % Generate offsprings and assign solutions to referece vectors
        for t = 1 : T            
            ChildDec = OffspringGeneration(Global,PopDec);
            ChildDec = unique(ChildDec,'rows','stable');
            ChildObj = PredObjs(ChildDec,Global,Model);
            
            if ~Global.fitcons
                [~,ChildCon]  = feval(Global.problem,Global.M,ChildDec);
            else
                ChildCon = PredCons(ChildDec,Global,Model);
            end

            PopDec        = [PopDec;ChildDec];
            PopObj        = [PopObj;ChildObj];
            PopCon        = [PopCon;ChildCon];

            if isempty(PopCon)
                PopCV = zeros(size(PopObj,1),1);
            else
                PopCV = abs(sum(max(0,PopCon),2));
            end

            CVeps = mean(PopCV)*(sum(PopCV==0)/numel(PopCV));

            % Adapt referece vectors
            if sum(PopCV<=CVeps) > 1
                Zmin            = min(PopObj(PopCV<=CVeps,:),[],1);
                Zmax            = max(PopObj(PopCV<=CVeps,:),[],1);
                if all(Zmin ~= Zmax)
                    Vi(1:Global.N,:) = V0.*repmat(Zmax-Zmin,size(V0,1),1);
                end
            else
                Zmin  = min(PopObj,[],1);
                Vi(1:Global.N,:) = V0;
            end
            if Global.numRVset > 1
                Vn = -Vi;
            else
                Vn = [];
            end
            if T > 1 && t < T
                [index1,SE1] = EnvironmentalSelection(PopObj,PopCV,CVeps,Vi,Zmin,Zmax,'S');
                if all(~isnan(Zmax))
                    [index2,SE2] = EnvironmentalSelection(PopObj,PopCV,CVeps,Vn,Zmin,Zmax,'D');
                else
                    index2 = [];
                    SE2    = inf;
                end
                if SE1 <= SE2
                    index  = index1;
                else
                    index  = index2;
                end
                PopDec = PopDec(index,:);
                PopObj = PopObj(index,:);
                if ~isempty(PopCon)
                    PopCon = PopCon(index,:);
                    PopCV = abs(sum(max(0,PopCon),2));
                else
                    PopCV = zeros(size(PopObj,1),1);
                end
            end
        end

        [ChildSelDec,RestDec,CVeps] = InfillSelection(Global,A2,Model,PopDec,PopObj,PopCV,CVeps,mu,Vi,Vn,Zmin,Zmax);
        NumOffsprings = [NumOffsprings;size(ChildSelDec,1)];

        if ~isempty(ChildSelDec)
            if ~Global.fitcons
                [~,ChildSelCon] = feval(Global.problem,Global.M,ChildSelDec);
            else
                ChildSelCon = PredCons(ChildSelDec,Global,Model);
            end
        else
            ChildSelCon = [];
        end
        if ~isempty(ChildSelCon)
            ChildSelCV = abs(sum(max(0,ChildSelCon),2));
        else
            ChildSelCV = zeros(size(ChildSelDec,1),1);
        end
        if ~isempty(RestDec)
            if ~Global.fitcons
                [~,RestCon]  = feval(Global.problem,Global.M,RestDec);
            else
                RestCon = PredCons(RestDec,Global,Model);
            end
        else
            RestCon = [];
        end
         
        if ~isempty(ChildSelDec)
            ChildSelObj  = nan(size(ChildSelDec,1),numel(Model.f));
            if ~isempty(ChildSelDec(ChildSelCV<=CVeps,:))
                New     = EvaluateSolutions(ChildSelDec(ChildSelCV<=CVeps,:),Global);
                A2.decs = [A2.decs;New.decs];
                A2.objs = [A2.objs;New.objs];
                A2.cons = [A2.cons;New.cons];
                
                f_train  = A2.objs;
                for m = 1 : Global.M
                    dmodel          = TrainMultiSurrogate(Global,A2.decs,f_train(:,m),THETAobj(m,:));
                    Model.f{m}      = dmodel;
                    if strcmpi(dmodel.surro_type,'krig')
                        THETAobj(m,:) = dmodel.theta;
                    end
                end
                if Global.fitcons
                    g_train	= A2.cons;
                    for n = 1 : numel(Model.g)
                        dmodel          = TrainMultiSurrogate(Global,A2.decs,g_train(:,n),THETAcon(n,:));
                        Model.g{n}      = dmodel;
                        if strcmpi(dmodel.surro_type,'krig')
                            THETAcon(n,:) = dmodel.theta;
                        end
                    end
                end
                if sum(ChildSelCV<=CVeps) > 0
                    ChildSelObj(ChildSelCV<=CVeps,:) = New.objs;
                end
                if sum(ChildSelCV>CVeps) > 0
                    ChildSelObj(ChildSelCV>CVeps,:) = PredObjs(ChildSelDec(ChildSelCV>CVeps,:),Global,Model);
                end
            else
                ChildSelObj = PredObjs(ChildSelDec,Global,Model);
            end
        else
            ChildSelObj = [];
        end

        RestObj        = PredObjs(RestDec,Global,Model);
        [id1,id2]      = ismember(RestDec,A2.decs,'rows');
        RestObj(id1,:) = A2.objs((id2(id2~=0)),:);

        PopAllDec = [ChildSelDec;RestDec];
        PopAllObj = [ChildSelObj;RestObj];
        PopAllCon = [ChildSelCon;RestCon];
        
        if isempty(PopAllCon)
            PopAllCV  = zeros(size(PopAllObj,1),1);
        else
            PopAllCV  = abs(sum(max(0,PopAllCon),2));
        end
                
        if sum(PopAllCV<=CVeps) > 1
            Zmin = min(PopAllObj(PopAllCV<=CVeps,:),[],1);
            Zmax = max(PopAllObj(PopAllCV<=CVeps,:),[],1);
            if all(Zmin ~= Zmax)
                Vi(1:Global.N,:) = V0.*repmat(Zmax-Zmin,size(V0,1),1);
            end
        else
            Zmin  = min(PopAllObj,[],1);
            Vi(1:Global.N,:) = V0;
        end
        if Global.numRVset > 1
            Vn = -Vi;
        else
            Vn = [];
        end
        TempAllDec = [PopAllDec;A2.decs];
        TempAllObj = [PopAllObj;A2.objs];
        TempAllCon = [PopAllCon;A2.cons];
        
        if ~isempty(TempAllCon)
            TempAllCV  = abs(sum(max(0,TempAllCon),2));
        else
            TempAllCV  = zeros(size(TempAllObj,1),1);
        end
        
        [index1,SE1] = EnvironmentalSelection(TempAllObj,TempAllCV,CVeps,Vi,Zmin,Zmax,'S');
        if all(~isnan(Zmax))
            [index2,SE2] = EnvironmentalSelection(TempAllObj,TempAllCV,CVeps,Vn,Zmin,Zmax,'D');
        else
            index2 = [];
            SE2    = inf;
        end
        
        if SE1 <= SE2
            NextIDs = index1;
        else
            NextIDs = index2;
        end
        
        clear A1;
        A1  = ReEvaluateSolutions(TempAllDec(NextIDs,:),A2,Global,Model);
        Gen = Gen + 1;
    end
    done = round(((size(A2.decs,1)/Global.FEmax)*100)*100)/100;
    disp([Global.problem,':  ',num2str(Global.M),'obj. ','Run-',num2str(Global.run),' Gen-',num2str(Gen),' (',num2str(done),' % finished)' ]);
    runtime = toc;
    Population = A2;
    save(['HSMEA_',Global.problem,'_M',num2str(Global.M),'_',num2str(Global.run),'.mat'],'Population','Global','NumOffsprings','runtime')
end

function Y = PredObjs(X,Global,Model)
try
    loop = numel(Model.f);
catch
    loop = numel(Model(1).f);
end
Y = nan(size(X,1),loop);
for i = 1:loop
    try
        model = Model{i};
    catch
        model = Model.f{i};
    end
    modeltype = model.surro_type;
    pred_func = strcat(modeltype, '_predict');
    Y(:,i) = feval(pred_func,Normalize(X,[Global.lower(:),Global.upper(:)]),model);
end
end

function Y = PredCons(X,Global,Model)
try
    loop = numel(Model.g);
catch
    loop = numel(Model(1).g);
end
Y = nan(size(X,1),loop);
for i = 1:loop
    try
        model = Model{i};
    catch
        model = Model.g{i};
    end
    modeltype = model.surro_type;
    pred_func = strcat(modeltype, '_predict');
    Y(:,i) = feval(pred_func,Normalize(X,[Global.lower(:),Global.upper(:)]),model);
end
end