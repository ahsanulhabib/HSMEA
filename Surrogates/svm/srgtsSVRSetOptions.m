function srgtOPT = srgtsSVRSetOptions(P, T, SVR_Kernel, SVR_KernelOptions, SVR_C, SVR_Loss, SVR_Insensitivity)
%Function srgtsSVRSetOptions creates the SURROGATES Toolbox option
%structure for support vector regression models. This structure contains
%the following fields:
%
%* GENERAL PARAMETERS
%
%   SRGT - Identifier of the surrogate technique: 'SVR'.
%   P    - NPOINTS-by-NDV matrix, where NPOINTS is the number of points of
%          the sample and NDV is the number of design variables.
%          Default: Empty matrix.
%   T    - NPOINTS-by-1 vector of responses on the P matrix points.
%          Default: Empty vector.
%
%* SUPPORT VECTOR REGRESSION PARAMETERS
%
%   SVR_Kernel        - Kernel [ string | 'Linear' | 'Polynomial' |
%                       'GaussianRBF' | 'MultiLayerPerceptron' |
%                       'LinearSpline' | 'LinearBSpline' |
%                       'TrigonometricPolynomial' | 'ExponentialRBF' |
%                       'AnovaSpline-1' | 'AnovaSpline-2' | 'AnovaSpline-3' |
%                       'AnovaBSpline' ]. Default: 'GaussianRBF'
%   SVR_KernelOptions - Kernel option (it depends on the kernel function):
%                         * Linear                   : no options required.
%                         * Polynomial               : degree.
%                         * GaussianRBF              : sigma. Default: 2.
%                         * MultiLayerPerceptron     : [p1 p2], where p1 is
%                                                    scale and p2 is
%                                                    offset.
%                         * LinearSpline             : no options required.
%                         * LinearBSpline            : degree.
%                         * TrigonometricPolynomial  : degree.
%                         * ExponentialRBF           : sigma.
%                         * AnovaSpline-1            : no options required.
%                         * AnovaSpline-2            : no options required.
%                         * AnovaSpline-3            : no options required.
%                         * AnovaBSpline             : maximum order of
%                                                    terms.
%
%   SVR_C             - Regularization parameter [ scalar ]. Default: Inf.
%   SVR_Loss          - Loss function [ string | 'einsensitive' |
%                       'quadratic' ]. Default: 'einsensitive'.
%   SVR_Insensitivity - Insensitivity [ scalar ]. Default: 0.
%
%The SURROGATES Toolbox uses the SVM Toolbox by Gunn (1997) to execute the
%support vector regression algorithm.
%
%This is how you can use srgtsSVRSetOptions:
%
%     OPTIONS = srgtsSVRSetOptions: creates a structure with the empty
%     parameters.
%
%     OPTIONS = srgtsSVRSetOptions(P, T): Given the sampled data P (input
%     variables) and T (output variables), it creates a structure with
%     default parameters used for all not specified fields.
%
%     OPTIONS = srgtsSVRSetOptions(P, T, ...
%         SVR_Kernel, SVR_KernelOptions, SVR_C, SVR_Loss,
%         SVR_Insensitivity)): it creates a structure with each of the
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
%     options = srgtsSVRSetOptions(X, Y)
%
%     options =
%
%                  SRGT: 'SVR'
%                     P: [5x1 double]
%                     T: [5x1 double]
%            SVR_Kernel: 'GaussianRBF'
%     SVR_KernelOptions: 2
%                 SVR_C: Inf
%              SVR_Loss: 'einsensitive'
%     SVR_Insensitivity: 0
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
srgtOPT.SRGT = 'SVR';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options
switch nargin
    case 0
        srgtOPT.P = [];
        srgtOPT.T = [];
        
        srgtOPT.FIT_Fn     = [];
        
        srgtOPT.SVR_Kernel        = [];
        srgtOPT.SVR_KernelOptions = [];
        srgtOPT.SVR_C             = [];
        srgtOPT.SVR_Loss          = [];
        srgtOPT.SVR_Insensitivity = [];
        
    case 2
        srgtOPT.P = P;
        srgtOPT.T = T;
        
        srgtOPT.FIT_Fn     = @svrgunn_fit;
        
        srgtOPT.SVR_Kernel        = 'GaussianRBF';
        srgtOPT.SVR_KernelOptions = 2;
        srgtOPT.SVR_C             = Inf;
        srgtOPT.SVR_Loss          = 'quadratic';
        srgtOPT.SVR_Insensitivity = 0;
    otherwise
        srgtOPT.P = P;
        srgtOPT.T = T;
        
        srgtOPT.FIT_Fn     = @svrgunn_fit;
        
        srgtOPT.SVR_Kernel        = SVR_Kernel;
        srgtOPT.SVR_KernelOptions = SVR_KernelOptions;
        srgtOPT.SVR_C             = SVR_C;
        srgtOPT.SVR_Loss          = SVR_Loss;
        srgtOPT.SVR_Insensitivity = SVR_Insensitivity;
end

return
