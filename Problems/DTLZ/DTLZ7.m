function [f,g] = DTLZ7(objnum,x)
if nargin <= 1
	prob.nf = objnum;
	prob.ng = 0;
    prob.k = 10-objnum+1;
    prob.nx = prob.nf + prob.k-1;
    for i = 1:prob.nx
        prob.range(i,:) = [0,1];
    end
    f = prob;
    g = [];
else
    [f,g] = dtlz7_true(x,objnum);
end
return

function [f,g] = dtlz7_true(x,objnum)
M = objnum;
k = 10-objnum+1;
f(:,1:M-1) = x(:,1:M-1);
G = 1+9/k*sum(x(:,M:M+k-1),2);
H = M - sum(f(:,1:M-1).*(1 + sin(3*pi*f(:,1:M-1))),2)./(1+G);
f(:,M) = (1+G).*H;
g = [];
return