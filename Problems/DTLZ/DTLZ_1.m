function [f,g] = DTLZ_1(objnum,x)
if nargin <= 1
	prob.nf = objnum;
	prob.ng = 0;
    prob.k = 10-objnum+1;
	prob.nx = prob.nf+prob.k-1;
	for i = 1:prob.nx
		prob.range(i,:) = [0,1];
	end
	f = prob;
else
	[f,g] = dtlz_1_true(x,objnum);
end
return

function [f,g] = dtlz_1_true(x,objnum)
M = objnum;
k = 10-objnum+1;
temp = sum((x(:,M:M+k-1)-0.5).^2-cos(20*pi*(x(:,M:M+k-1)-0.5)),2);
G = 100*(k+temp);
f(:,1) = 0.5*prod(x(:,1:objnum-1),2).*(1+G);
for i = 2:M-1
    f(:,i) = 0.5*prod(x(:,1:M-i),2).*(1-x(:,M-i+1)).*(1+G);
end
f(:,M) = 0.5*(1-x(:,1)).*(1+G);
g = [];
f = -f;
return