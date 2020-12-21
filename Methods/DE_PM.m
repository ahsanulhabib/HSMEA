function Offspring = DE_PM(Global,ParentDec)
% Differential Evolution followed by Polynomial Mutation
% Modified by Ahsanul Habib
%------------------------------- Copyright --------------------------------
% Copyright (c) 2016-2017 BIMK Group.
%------------------------------- Reference --------------------------------
% "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    [N,D]     = size(ParentDec);
    [CR,F,proM,disM] = deal(Global.CR,Global.F,1/D,Global.mutation_poly);

    %% Differental evolution
    Parent1Dec   = ParentDec(1:N/3,:);
    Parent2Dec   = ParentDec(N/3+1:N/3*2,:);
    Parent3Dec   = ParentDec(N/3*2+1:end,:);
    OffspringDec = Parent1Dec;
    Site = rand(N/3,D) < CR;
    OffspringDec(Site) = OffspringDec(Site) + F*(Parent2Dec(Site)-Parent3Dec(Site));

    %% Polynomial mutation
    Site  = rand(N/3,D) < proM/D;
    mu    = rand(N/3,D);
    temp1 = Site & mu<=0.5;
    Lower = repmat(Global.lower,N/3,1);
    Upper = repmat(Global.upper,N/3,1);
    OffspringDec(temp1) = OffspringDec(temp1)+(Upper(temp1)-Lower(temp1)).*((2.*mu(temp1)+(1-2.*mu(temp1)).*...
                         (1-(OffspringDec(temp1)-Lower(temp1))./(Upper(temp1)-Lower(temp1))).^(disM+1)).^(1/(disM+1))-1);
    temp2  = Site & mu>0.5; 
    OffspringDec(temp2) = OffspringDec(temp2)+(Upper(temp2)-Lower(temp2)).*(1-(2.*(1-mu(temp2))+2.*(mu(temp2)-0.5).*...
                         (1-(Upper(temp2)-OffspringDec(temp2))./(Upper(temp2)-Lower(temp2))).^(disM+1)).^(1/(disM+1)));
                     
    temp3 = OffspringDec < Lower;
    temp4 = OffspringDec > Upper;
    
    OffspringDec(temp3) = Lower(temp3);
    OffspringDec(temp4) = Upper(temp4);
    
    Offspring           = OffspringDec(randi(size(OffspringDec,1)),:);
end