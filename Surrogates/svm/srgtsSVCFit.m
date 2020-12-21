function srgtSRGT = srgtsSVCFit(srgtOPT)
%Function srgtsSVCFit fits the specified support vector classification model
%using the SVM Toolbox of Gunn (1997).
% 
%    srgtSRGT = srgtsSVCFit(srgtOPT)
% 
%srgtSRGT is the surrogate structure that contains the following fields:
%* P                 : experimental design matrix.
%* SVR_Kernel        : Kernel function.
%* SVR_KernelOptions : options for the kernel function.
%* SVR_NbSV          : number of support vectors.
%* SVR_DiffLagMult   : difference of Lagrange multipliers.
%* SVR_Bias          : bias term.
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
%     options = srgtsSVRSetOptions(X, Y);
% 
%     surrogate = srgtsSVCFit(options)
% 
%     surrogate = 
% 
%                     P: [5x1 double]
%            SVC_Kernel: 'GaussianRBF'
%     SVC_KernelOptions: 2
%              SVC_NbSV: 5
%       SVC_DiffLagMult: [5x1 double]
%              SVC_Bias: 0
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

srgtSRGT.P                 = srgtOPT.P;
srgtSRGT.T                 = srgtOPT.T;
srgtSRGT.SVC_Kernel        = srgtOPT.SVC_Kernel;
srgtSRGT.SVC_KernelOptions = srgtOPT.SVC_KernelOptions;

[srgtSRGT.SVC_NbSV, srgtSRGT.SVC_LagMult, srgtSRGT.SVC_Bias] = svcgunn_fit(...
    srgtOPT.P, ...
    srgtOPT.T, ...
    srgtOPT.SVC_Kernel, ...
    srgtOPT.SVC_KernelOptions, ...
    srgtOPT.SVC_C);

return
