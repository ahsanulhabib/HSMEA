function Offspring = SBX_PM(Global,ParentDec)
% Simulated Binary Crossover followed by Polynomial Mutation
% Modified by Ahsanul Habib
%------------------------------- Copyright --------------------------------
% Copyright (c) 2016-2017 BIMK Group.
%------------------------------- Reference --------------------------------
% "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    [N,D]     = size(ParentDec);
    [proC,disC,proM,disM] = deal(Global.crossover_pr,Global.crossover_sbx,1/D,Global.mutation_poly);
   
    %% Simulated binary crossover
    Parent1Dec = ParentDec(1:N/2,:);
    Parent2Dec = ParentDec(N/2+1:end,:);
    beta = zeros(N/2,D);
    mu   = rand(N/2,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N/2,D);
    beta(rand(N/2,D)<0.5) = 1;
    beta(repmat(rand(N/2,1)>proC,1,D)) = 1;
    OffspringDec = [(Parent1Dec+Parent2Dec)/2+beta.*(Parent1Dec-Parent2Dec)/2
                    (Parent1Dec+Parent2Dec)/2-beta.*(Parent1Dec-Parent2Dec)/2];

    %% Polynomial mutation
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp1 = Site & mu<=0.5;
    Lower = repmat(Global.lower,N,1);
    Upper = repmat(Global.upper,N,1);
    OffspringDec(temp1) = OffspringDec(temp1)+(Upper(temp1)-Lower(temp1)).*((2.*mu(temp1)+(1-2.*mu(temp1)).*...
                         (1-(OffspringDec(temp1)-Lower(temp1))./(Upper(temp1)-Lower(temp1))).^(disM+1)).^(1/(disM+1))-1);
    temp2 = Site & mu>0.5; 
    OffspringDec(temp2) = OffspringDec(temp2)+(Upper(temp2)-Lower(temp2)).*(1-(2.*(1-mu(temp2))+2.*(mu(temp2)-0.5).*...
                         (1-(Upper(temp2)-OffspringDec(temp2))./(Upper(temp2)-Lower(temp2))).^(disM+1)).^(1/(disM+1)));
                     
    temp3 = OffspringDec < Lower;
    temp4 = OffspringDec > Upper;
    
    OffspringDec(temp3) = Lower(temp3);
    OffspringDec(temp4) = Upper(temp4);
    
    Offspring           = OffspringDec(randi(size(OffspringDec,1)),:);
end