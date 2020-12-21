function [f,g] = WFG9(objnum,x)
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
    [f,g] = wfg9_true(x,objnum);
end
return

function [f,g] = wfg9_true(x,objnum)
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
Y = (fliplr(cumsum(fliplr(z01),2))-z01)./repmat(K+L-1:-1:0,N,1);
t1(:,1:K+L-1) = z01(:,1:K+L-1).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,1:K+L-1)).*abs(floor(0.5-Y(:,1:K+L-1))+0.98/49.98)));
t1(:,end)     = z01(:,end);

t2 = zeros(N,K+L);
t2(:,1:K)     = s_decept(t1(:,1:K),0.35,0.001,0.05);
t2(:,K+1:end) = s_multi(t1(:,K+1:end),30,95,0.35);

t3 = zeros(N,M);
for i = 1 : M-1
    t3(:,i) = r_nonsep(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),K/(M-1));
end

SUM = zeros(N,1);
for i = K+1 : K+L-1
    for j = i+1 : K+L
        SUM = SUM + abs(t2(:,i)-t2(:,j));
    end
end

t3(:,M) = (sum(t2(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));

tempf = zeros(N,M);
for i = 1 : M-1
    tempf(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
end
tempf(:,M) = t3(:,M);

h = concave(tempf);
f = repmat(D*tempf(:,M),1,M) + repmat(S,N,1).*h;

g = [];
return

function Output = s_decept(y,A,B,C)
Output = 1+(abs(y-A)-B).*(floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
return

function Output = s_multi(y,A,B,C)
Output = (1+cos((4*A+2)*pi*(0.5-abs(y-C)/2./(floor(C-y)+C)))+4*B*(abs(y-C)/2./(floor(C-y)+C)).^2)/(B+2);
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