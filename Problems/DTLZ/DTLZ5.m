function [f,g] = DTLZ5(objnum,x)
if nargin <= 1
	prob.nf = objnum;
	prob.ng = 0;
    prob.k = 10-objnum+1;
	prob.nx = prob.nf+prob.k-1;
	for i = 1:prob.nx
		prob.range(i,:) = [0,1];
	end
	f = prob;
    g = [];
else
	[f,g] = dtlz5_true(objnum,x);
end
return

function [f,g] = dtlz5_true(objnum,x)
M = objnum;
k = 10-objnum+1;
temp = x(:,M:M+k-1)-0.5;
G = sum(temp.^2,2);
temp1(:,1) = x(:,1)*pi/2;
for i = 1:numel(G)
    temp1(i,2:M-1) = pi*(1+2*G(i)*x(i,2:M-1))/(4*(1+G(i)));
end
f(:,1) = (1+G).*prod(cos(temp1),2);
for j = 2:M-1
    f(:,j)= (1+G).*prod(cos(temp1(:,1:M-j)),2).*sin(temp1(:,M-j+1));
end
f(:,M) = (1+G).*sin(temp1(:,1));
g = [];
return