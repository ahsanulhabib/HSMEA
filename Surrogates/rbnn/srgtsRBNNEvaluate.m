function yhat = srgtsRBNNEvaluate(x, srgtSRGT)
%Function srgtsRBNNEvaluate is used to predict the response of a radial
%basis neural network model. Thus, for example:
%
%     YHAT = srgtsRBNNEvaluate(X, SURROGATE): returns the response YHAT
%     predicted by the radial basis neural network model SURROGATE at all X
%     sites. X can be either a single row vector (single point, with each
%     column representing a variable) or a matrix (each row represents a
%     point). YHAT is an NBPOINTS-by-1 vector.
%
%Example:
%     % basic information about the problem
%     myFN = @cos;  % this could be any user-defined function
%      = [0;     % lower bound
%                    2*pi]; % upper bound
%
%     % create DOE
%     npoints = 5;
%     X = linspace((1), (2), npoints)';
%
%     % evaluate analysis function at X points
%     Y = feval(myFN, X);
%
%     % fit surrogate models
%     options = srgtsRBNNSetOptions(X, Y);
% 
%     surrogate = srgtsRBNNFit(options);
%
%     % create test points
%     Xtest = linspace((1), (2), 100)';
%
%     % evaluate surrogate at Xtest
%     Yhat = srgtsRBNNEvaluate(Xtest, surrogate);
%
%     plot(X, Y, 'o', ...
%          Xtest, Yhat)

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

yhat = sim(srgtSRGT.RBNN_NEWRBModel, x.');
yhat = yhat.';

return
