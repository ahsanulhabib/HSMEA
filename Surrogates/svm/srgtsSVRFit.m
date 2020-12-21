function [srgtSRGT,srgtSTT] = srgtsSVRFit(srgtOPT)
%Function srgtsSVRFit fits the specified support vector regression model
%using the SVM Toolbox of Gunn (1997).
%
%    [srgtSRGT srgtSTT] = srgtsSVRFit(srgtOPT)
%
%srgtSRGT is the surrogate structure that contains the following fields:
%* P                 : experimental design matrix.
%* SVR_Kernel        : Kernel function.
%* SVR_KernelOptions : options for the kernel function.
%* SVR_NbSV          : number of support vectors.
%* SVR_DiffLagMult   : difference of Lagrange multipliers.
%* SVR_Bias          : bias term.
%
%srgtSTT is the state structure that contains the following fields:
%* FIT_Fn       : function handle of the fitting function.
%* FIT_LossFn   : function handle of the loss function.
%* FIT_ObjVal   : value of the loss function (after fitting).
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
%     surrogate = srgtsSVRFit(options)
%
%     surrogate =
%
%                     P: [5x1 double]
%            SVR_Kernel: 'GaussianRBF'
%     SVR_KernelOptions: 2
%              SVR_NbSV: 5
%       SVR_DiffLagMult: [5x1 double]
%              SVR_Bias: 0
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
srgtSRGT.SVR_Kernel        = srgtOPT.SVR_Kernel;
srgtSRGT.SVR_KernelOptions = srgtOPT.SVR_KernelOptions;

srgtSTT = srgtsFitCreateState(srgtOPT);

[srgtSRGT.SVR_NbSV, srgtSRGT.SVR_DiffLagMult, srgtSRGT.SVR_Bias] = svrgunn_fit(...
                                                                                srgtOPT.P, ...
                                                                                srgtOPT.T, ...
                                                                                srgtOPT.SVR_Kernel, ...
                                                                                srgtOPT.SVR_KernelOptions, ...
                                                                                srgtOPT.SVR_C, ...
                                                                                srgtOPT.SVR_Loss, ...
                                                                                srgtOPT.SVR_Insensitivity);
return
