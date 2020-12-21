function Score = SEnergy(PopObj)
[N,M] = size(PopObj);
Score = 0;
for i = 1:N-1
    j = i+1:N;
    Distances = pdist2(PopObj(i,:),(PopObj(j,:)));
    Distances(Distances==0) = 1e-6;
    Score = Score + sum((Distances).^-(M-1));
end
end