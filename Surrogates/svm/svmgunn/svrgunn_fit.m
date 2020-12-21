function [nsv, beta, bias] = svrgunn_fit(X, Y, ...
    SVR_Kernel, SVR_KernelOptions, SVR_C, SVR_Loss, SVR_Insensitivity)
%SVR Support Vector Regression
%
%  Usage: [nsv beta bias] = svr(X,Y,SVR_Kernel,SVR_KernelOptions,SVR_C,SVR_Loss,SVR_Insensitivity)
%
%  Parameters: X          - Training inputs
%              Y          - Training targets
%              SVR_Kernel        - kernel function
%              SVR_KernelOptions - kernel function
%              SVR_C          - upper bound (non-separable case)
%              SVR_Loss       - SVR_Loss function
%              SVR_Insensitivity          - insensitivity
%              nsv        - number of support vectors
%              beta       - Difference of Lagrange Multipliers
%              bias       - bias term
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tolerance for support vector detection
% epsilon = svtol(SVR_C);
if SVR_C == Inf
    epsilon = 1e-5;
else
    epsilon = SVR_C*1e-6;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the kernel matrix
n = size(X,1);
H = zeros(n,n);
for i=1:n
    for j=1:n
        H(i,j) = svmgunn_kernelmatrix(SVR_Kernel,SVR_KernelOptions,X(i,:),X(j,:));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the parameters for the optimization problem
switch lower(SVR_Loss)
    case 'einsensitive',
        Hb = [H -H; -H H];
        c  = [ (SVR_Insensitivity*ones(n,1) - Y) ;
            (SVR_Insensitivity*ones(n,1) + Y) ];

        vlb = zeros(2*n,1);  % set the bounds: alphas >= 0
        vub = SVR_C*ones(2*n,1); %                 alphas <= SVR_C
        x0  = zeros(2*n,1);  % the starting point is [0 0 0 ... 0]
        
        neqcstr = svmgunn_nobias(SVR_Kernel); % set the number of equality constraints (1 or 0)
        if neqcstr
            A = [ones(1,n) -ones(1,n)];
            b = 0; % set the constraint Ax = b
        else
            A = [];
            b = [];
        end
        
    case 'quadratic',
        Hb = H + eye(n)/(2*SVR_C);
        c  = -Y;
        
        vlb = -1e30*ones(n,1);
        vub = 1e30*ones(n,1);
        x0  = zeros(n,1); % the starting point is [0 0 0   0]
        
        neqcstr = svmgunn_nobias(SVR_Kernel); % set the number of equality constraints (1 or 0)
        if neqcstr
            A = ones(1,n);
            b = 0; % set the constraint Ax = b
        else
            A = [];
            b = [];
        end
        
    otherwise
        disp('Error: Unknown Loss Function\n');
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add small amount of zero order regularisation to avoid problems when Hessian
% is badly conditioned. Rank is always less than or equal to n. Note that adding
% to much reg will peturb solution

Hb = Hb + 1e-10*eye(size(Hb));

% solve the optimisation problem
if ~strcmpi(SVR_Loss,'quadratic')
    alpha = svmgunn_qp(Hb, c, A, b, vlb, vub, x0, neqcstr);
else
    options = optimoptions('quadprog','Algorithm','interior-point-convex','display','off');
    alpha = quadprog(Hb, c, A, b, [], [], vlb, vub, x0, options);
end

switch lower(SVR_Loss)
    case 'einsensitive',
        beta =  alpha(1:n) - alpha(n+1:2*n);
        
    case 'quadratic',
        beta = alpha;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the number of support vectors
svi = find( abs(beta) > epsilon );
nsv = length( svi );

% implicit bias, b0
bias = 0;

% Explicit bias, b0
if svmgunn_nobias(SVR_Kernel) ~= 0
    
    switch lower(SVR_Loss)
        case 'einsensitive',
            % find bias from average of support vectors with interpolation
            % error. SVs with interpolation error SVR_Insensitivity have
            % alphas: 0 < alpha < SVR_C
            svii = find( abs(beta) > epsilon & abs(beta) < (SVR_C - epsilon));
            
            if ~isempty(svii)    
                bias = (1/length(svii))*sum(Y(svii) - SVR_Insensitivity*sign(beta(svii)) - H(svii,svi)*beta(svi));
            else
                fprintf('No support vectors with interpolation error SVR_Insensitivity - cannot compute bias.\n');
                bias = (max(Y)+min(Y))/2;
            end
            
        case 'quadratic',
            bias = mean(Y - H*beta);
    end
    
end

return
