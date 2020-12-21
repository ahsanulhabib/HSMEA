function [nsv, alpha, b0] = svcgunn_fit(X, Y, SVC_Kernel, SVC_KernelOptions, SVC_C)
%SVC Support Vector Classification
%
%  Usage: [nsv alpha bias] = svcgunn_fit(X,Y,SVC_Kernel,SVC_C)
%
%  Parameters: X      - Training inputs
%              Y      - Training targets
%              SVR_Kernel        - kernel function
%              SVR_KernelOptions - kernel function
%              SVR_C          - upper bound (non-separable case)
%              nsv    - number of support vectors
%              alpha  - Lagrange Multipliers
%              b0     - bias term
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

n = size(X,1);

% tolerance for Support Vector Detection
if SVC_C == Inf
    epsilon = 1e-5;
else
    epsilon = SVC_C*1e-6;
end

% Construct the Kernel matrix
H = zeros(n,n);
for i=1:n
    for j=1:n
        H(i,j) = Y(i)*Y(j)*svmgunn_kernelmatrix(SVC_Kernel, SVC_KernelOptions, X(i,:), X(j,:));
    end
end
c = -ones(n,1);

% Add small amount of zero order regularisation to
% avoid problems when Hessian is badly conditioned.
H = H+1e-10*eye(size(H));

% Set up the parameters for the Optimisation problem

vlb = zeros(n,1);      % Set the bounds: alphas >= 0
vub = SVC_C*ones(n,1);     %                 alphas <= SVC_C
x0  = zeros(n,1);       % The starting point is [0 0 0   0]
neqcstr = svmgunn_nobias(SVC_Kernel); % Set the number of equality constraints (1 or 0)
if neqcstr
    A = Y';
    b = 0;     % Set the constraint Ax = b
else
    A = [];
    b = [];
end

% Solve the Optimisation Problem
alpha = svmgunn_qp(H, c, A, b, vlb, vub, x0, neqcstr);

% Compute the number of Support Vectors
svi = find( alpha > epsilon);
nsv = length(svi);

% Implicit bias, b0
b0 = 0;

% Explicit bias, b0
if svmgunn_nobias(SVC_Kernel) ~= 0
    % find b0 from average of support vectors on margin
    % SVs on margin have alphas: 0 < alpha < SVC_C
    svii = find( alpha > epsilon & alpha < (SVC_C - epsilon));
    if ~isempty(svii)
        b0 =  (1/length(svii))*sum(Y(svii) - H(svii,svi)*alpha(svi).*Y(svii));
    else
        fprintf('No support vectors on margin - cannot compute bias.\n');
    end
end

return
