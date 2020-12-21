function [f,g] = DTLZ6(objnum,x)
if nargin <= 1
	prob.nf = objnum;
	prob.ng = 0;
    prob.k = 10-objnum+1;
    prob.nx =  prob.nf+prob.k-1;
    for i = 1:prob.nx
        prob.range(i,:) = [0,1];
    end
    f = prob;
    g = [];
else
    [f,g] = dtlz6_true(x,objnum);
end
return

function [f,g] = dtlz6_true(x,objnum)
M = objnum;
k = 10-objnum+1;
temp = x(:,M:M+k-1);
G = sum(temp.^0.1,2);
theta(:,1) = x(:,1)*pi/2;
for i = 1:numel(G)
    theta(i,2:M-1) = pi*(1+2*G(i)*x(i,2:M-1))/(4*(1+G(i)));
end
f(:,1) = (1+G).*prod(cos(theta),2);
for j = 2:M-1
    f(:,j) = (1+G).*prod(cos(theta(:,1:M-j)),2).*sin(theta(:,M-j+1));
end
f(:,M) = (1+G).*sin(theta(:,1));
g = []; 
return