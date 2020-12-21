function Generate_WFGminus_PFs
parfor i = 1:9
    Objs = [3,4,6,8,10];

    for j = 1:numel(Objs)
        M = Objs(j);
        h = load(['WFG',num2str(i),'.',num2str(M),'D.pf',]);
        Range = minmax(h')';
        PFt = (h - repmat(Range(1,:),size(h,1),1)) ./ repmat((Range(2,:)-Range(1,:)),size(h,1),1);
        PF  = repmat(-2*(1:M),size(h,1),1) + repmat((-1*ones(1,M)+2*(1:M)),size(h,1),1).*(1-PFt);
        savepf(['WFG_',num2str(i)],Objs(j),PF);
    end
end
end

function savepf(f_name,M,PF)
save([f_name,'.',num2str(M),'D.pf'],'PF','-ascii');
end 