function srgtOPT = srgtsSVCSetOptions(P, T, SVC_Kernel, SVC_KernelOptions, SVC_C)
%Function srgtsSVCSetOptions creates the SURROGATES Toolbox option
%structure for support vector classification models. This structure
%contains the following fields:
%
%* GENERAL PARAMETERS
%
%   SRGT - Identifier of the surrogate technique: 'SVC'.
%   P    - NPOINTS-by-NDV matrix, where NPOINTS is the number of points of
%          the sample and NDV is the number of design variables.
%          Default: Empty matrix.
%   T    - NPOINTS-by-1 vector of responses on the P matrix points.
%          Default: Empty vector.
%
%* SUPPORT VECTOR REGRESSION PARAMETERS
%
%   SVC_Kernel        - Kernel [ string | 'Linear' | 'Polynomial' |
%                       'GaussianRBF' | 'MultiLayerPerceptron' |
%                       'LinearSpline' | 'LinearBSpline' |
%                       'TrigonometricPolynomial' | 'ExponentialRBF' |
%                       'AnovaSpline-1' | 'AnovaSpline-2' | 'AnovaSpline-3' |
%                       'AnovaBSpline' ]. Default: 'GaussianRBF'
%   SVC_KernelOptions - Kernel option (it depends on the kernel function):
%                         * Linear                  : no options required.
%                         * Polynomial              : degree.
%                         * GaussianRBF             : sigma. Default: 2.
%                         * MultiLayerPerceptron    : [p1 p2], where p1 is
%                                                   scale and p2 is offset.
%                         * LinearSpline            : no options required.
%                         * LinearBSpline           : degree.
%                         * TrigonometricPolynomial : degree.
%                         * ExponentialRBF           : sigma.
%                         * AnovaSpline-1            : no options required.
%                         * AnovaSpline-2            : no options required.
%                         * AnovaSpline-3            : no options required.
%                         * AnovaBSpline             : maximum order of
%                                                    terms.
%
%   SVC_C             - Regularization parameter [ scalar ]. Default: Inf.
%
%The SURROGATES Toolbox uses the SVM Toolbox by Gunn (1997) to execute the
%support vector classification algorithm.
%
%This is how you can use srgtsSVCSetOptions:
% 
%     OPTIONS = srgtsSVCSetOptions: creates a structure with the empty
%     parameters.
%
%     OPTIONS = srgtsSVCSetOptions(P, T): Given the sampled data P (input
%     variables) and T (output variables), it creates a structure with
%     default parameters used for all not specified fields.
%
%     OPTIONS = srgtsSVCSetOptions(P, T, ...
%         SVC_Kernel, SVC_KernelOptions, SVC_C, SVC_Loss,
%         SVC_Insensitivity)): it creates a structure with each of the
%         specified fields.
%
%Example:
%     % basic information about the problem
%     myFN = @cos;  % this could be any user-defined function
%     designspace = [0;     % lower bound
%                    2*pi]; % upper bound
%
%     % create DOE
%     npoints = 5;
%     X = linspace(designspace(1), designspace(2), npoints)';
%
%     % evaluate analysis function at X points
%     Y = feval(myFN, X);
%
%     % fit surrogate models
%     options = srgtsSVCSetOptions(X, Y)
%
%     options =
%
%                  SRGT: 'SVC'
%                     P: [5x1 double]
%                     T: [5x1 double]
%            SVC_Kernel: 'GaussianRBF'
%     SVC_KernelOptions: 2
%                 SVC_C: Inf
%              SVC_Loss: 'einsensitive'
%     SVC_Insensitivity: 0
%
%REFERENCES:
%
%Gunn SR, Support Vector Machines for Classification and Regression,
%Technical Report, Image Speech and Intelligent Systems Research Group,
%University of Southampton, UK, 1997.
%Available at: http://www.isis.ecs.soton.ac.uk/resources/svminfo/.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
srgtOPT.SRGT = 'SVC';

srgtOPT.P = P;
srgtOPT.T = T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% svr options
if nargin == 2
    srgtOPT.SVC_Kernel        = 'GaussianRBF';
    srgtOPT.SVC_KernelOptions = 2;
    srgtOPT.SVC_C             = Inf;
else
    srgtOPT.SVC_Kernel        = SVC_Kernel;
    srgtOPT.SVC_KernelOptions = SVC_KernelOptions;
    srgtOPT.SVC_C             = SVC_C;
end

return
