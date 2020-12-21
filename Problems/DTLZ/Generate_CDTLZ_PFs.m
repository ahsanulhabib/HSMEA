function Generate_CDTLZ_PFs
for i = 1:3
    if i == 1 
        f_name = ['C1_DTLZ',num2str(i)];
    elseif i == 2
        f_name = ['C2_DTLZ',num2str(i)];
    elseif i == 3
        f_name = ['C3_DTLZ',num2str(i+1)];
    end
    Objs = [3,4,6,8,10];
    Pts = [5050,10660,33649,50388,92378];

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