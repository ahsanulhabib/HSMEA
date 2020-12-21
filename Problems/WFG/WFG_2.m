function [f,g] = WFG_2(objnum,x)
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
	[f,g] = wfg_2_true(x,objnum);
end
return

function [f,g] = wfg_2_true(x,objnum)
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

t2 = zeros(N,K+L/2);
t2(:,1:K) = t1(:,1:K);
t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;

t3 = zeros(N,M);
for i = 1 : M-1
    t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
end
t3(:,M) = r_sum(t2(:,K+1:K+L/2),ones(1,L/2));

tempf = zeros(N,M);
for i = 1 : M-1
    tempf(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
end
tempf(:,M) = t3(:,M);

h      = convex(tempf);
h(:,M) = disc(tempf);
f = -(repmat(D*tempf(:,M),1,M) + repmat(S,N,1).*h);

g = [];
return

function Output = s_linear(y,A)
Output = abs(y-A)./abs(floor(A-y)+A);
return

function Output = r_sum(y,w)
Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
return

function Output = convex(x)
Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
return

function Output = disc(x)
Output = 1-x(:,1).*(cos(5*pi*x(:,1))).^2;
return