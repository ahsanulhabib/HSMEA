function Yhat = svrgunn_evaluate(Ptest, Ptraining, SVR_Kernel, SVR_KernelOptions, SVR_DiffLagMult, SVR_Bias)
%SVRMSE Calculate SVR Output
%
%  Usage: Yhat = svrgunn_evaluate(Ptraining,Ptest,SVR_Kernel,SVR_KernelOptions,SVR_DiffLagMult,SVR_Bias)
%
%  Parameters: Ptraining         - Training inputs
%              Ptest             - Test inputs
%              SVR_Kernel        - kernel function
%              SVR_KernelOptions - kernel options
%              SVR_DiffLagMult   - Difference of Lagrange Multipliers
%              SVR_Bias          - SVR_Bias
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

n = size(Ptraining,1);
m = size(Ptest,1);
H = zeros(m,n);
for i=1:m
    for j=1:n
        H(i,j) = svmgunn_kernelmatrix(SVR_Kernel, SVR_KernelOptions, Ptest(i,:), Ptraining(j,:));
    end
end
Yhat = (H*SVR_DiffLagMult + SVR_Bias);
return
