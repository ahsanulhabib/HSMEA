function Yhat = svcgunn_evaluate(Ptest, Ptraining, Ytraining, SVC_Kernel, SVC_KernelOptions, SVC_LagMult, SVC_Bias, actfunc)
%SVCOUTPUT Calculate SVC Output
%
%  Usage: Yhat = svcgunn_evaluate(Ptraining, Ytraining,Ptest,SVC_Kernel,SVC_LagMult,SVC_Bias,actfunc)
%
%  Parameters: Ptraining   - Training inputs
%              Ytraining   - Training targets
%              Ptest   - Test inputs
%              SVC_Kernel        - kernel function
%              SVC_KernelOptions - kernel options
%              SVC_LagMult   - Lagrange Multipliers
%              SVC_Bias   - SVC_Bias
%              actfunc- activation function (0(default) hard | 1 soft)
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

if (nargin == 7)
    actfunc = 0;
end
n = size(Ptraining,1);
m = size(Ptest,1);
H = zeros(m,n);
for i=1:m
    for j=1:n
        H(i,j) = Ytraining(j)*svmgunn_kernelmatrix(SVC_Kernel, SVC_KernelOptions, Ptest(i,:), Ptraining(j,:));
    end
end
if (actfunc)
    Yhat = softmargin(H*SVC_LagMult + SVC_Bias);
else
    Yhat = sign(H*SVC_LagMult + SVC_Bias);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friend functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = softmargin(x)
%SOFTMARGIN Support Vector Classification Softmargin
%
%  Usage: y = softmargin(x)
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

y = x;
y(find(x < -1)) = -1;
y(find(x > 1)) = 1;

return
