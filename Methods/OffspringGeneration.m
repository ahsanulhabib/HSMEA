function ChildDec = OffspringGeneration(Global,PopDec)
ChildDec = [];
while size(ChildDec) < Global.N
    if rand < 0.5
        IDsbx = randi(size(PopDec,1),[2,1]);
        ChildDec = [ChildDec;SBX_PM(Global,PopDec(IDsbx,:))];
    else
        IDde = randi(size(PopDec,1),[3,1]);
        ChildDec = [ChildDec;DE_PM(Global,PopDec(IDde,:))];
    end
end
end