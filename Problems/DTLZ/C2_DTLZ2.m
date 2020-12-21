function [f,g] = C2_DTLZ2(objnum,x)
if nargin <= 1
	prob.nf = objnum;
	prob.ng = 1;
    prob.k = 10-objnum+1;    
	prob.nx = prob.nf+prob.k-1;
	for i = 1:prob.nx
		prob.range(i,:) =  [0,1];
	end
	f = prob;
else
	[f,g] = c2dtlz2_true(x,objnum);
end
return

function [f,g] = c2dtlz2_true(x,objnum)
M = objnum;
k = 10-objnum+1;
temp = x(:,M:M+k-1)-0.5;
G = sum(temp.^2,2);
temp1 = x(:,1:M-1)*pi/2;
f(:,1) = (1+G).*prod(cos(temp1),2);
for j = 2:M-1
    f(:,j)= (1+G).*prod(cos(temp1(:,1:M-j)),2).*sin(temp1(:,M-j+1));
end
f(:,M) = (1+G).*sin(x(:,1)*pi/2);
if M == 3
    r = 0.4;
else
    r = 0.5;
end
temp2 = [];
for i = 1:M
    id = [];
    id = setdiff(1:M,i);
    temp2(:,i) = sum(f(:,id).^2,2) - repmat(r.^2,size(f,1),1);
end
g = min([min((f - 1).^2 + temp2,[],2) sum((f - 1/sqrt(M)).^2,2) - repmat(r.^2,size(f,1),1)],[],2);
 return