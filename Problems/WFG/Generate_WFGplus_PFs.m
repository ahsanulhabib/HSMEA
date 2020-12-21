function Generate_WFGplus_PFs
for i = 1:9
    Objs = [3,4,6,8,10];
    Pts = [5050,10660,33649,50388,92378];

    f_name = ['WFG',num2str(i)];
    for j = 1:numel(Objs)
        PF = calculatePF(f_name,Objs(j),Pts(j));
        savepf(f_name,Objs(j),PF);
    end
end
end

function PF = calculatePF(f_name,M,Pts)
PF = feval(f_name,M,Pts);
end

function savepf(f_name,M,PF)
save([f_name,'.',num2str(M),'D.pf'],'PF','-ascii');
end