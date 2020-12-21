function [f,g] = WFG6(objnum,x)
if nargin <= 1
	prob.nf = objnum;
	prob.ng = 0;
    if prob.nf == 4
        prob.k = 6;
    else
        prob.k = objnum-1;
    end
    if prob.nf == 6
        prob.l = 9-prob.k;
    elseif prob.nf == 8
        prob.l = 9-prob.k;
    elseif prob.nf == 10
        prob.l = 11-prob.k;
    else
        prob.l = 10-prob.k;
    end
	prob.nx = prob.k+prob.l;
	for i = 1:prob.nx
		prob.range(i,:) = [0,2*i];
	end
	f = prob;
else
	[f,g] = wfg6_true(x,objnum);
end
return

function [f,g] = wfg6_true(x,objnum)
M = objnum;
if M == 4
    K = 6;
else
    K = objnum-1;
end
L = size(x,2) - K;
D = 1;
S = 2 : 2 : 2*M;
A = ones(1,M-1);
N = size(x,1);

z01 = x./repmat(2:2:size(x,2)*2,N,1);

t1 = zeros(N,K+L);
t1(:,1:K)     = z01(:,1:K);
t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);

t2 = zeros(N,M);
for i = 1 : M-1
    t2(:,i) = r_nonsep(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),K/(M-1));
end

SUM = zeros(N,1);
for i = K+1 : K+L-1
    for j = i+1 : K+L
        SUM = SUM + abs(t1(:,i)-t1(:,j));
    end
end

t2(:,M) = (sum(t1(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));

tempf = zeros(N,M);
for i = 1 : M-1
    tempf(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
end
tempf(:,M) = t2(:,M);

h = concave(tempf);
f = repmat(D*tempf(:,M),1,M) + repmat(S,N,1).*h;

g = [];
return

function Output = s_linear(y,A)
Output = abs(y-A)./abs(floor(A-y)+A);
return

function Output = r_nonsep(y,A)
    Output = zeros(size(y,1),1);
    for j = 1 : size(y,2)
        Temp = zeros(size(y,1),1);
        for k = 0 : A-2
            Temp = Temp+abs(y(:,j)-y(:,1+mod(j+k,size(y,2))));
        end
        Output = Output+y(:,j)+Temp;
    end
    Output = Output./(size(y,2)/A)/ceil(A/2)/(1+2*A-2*ceil(A/2));
return

function Output = concave(x)
Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
return