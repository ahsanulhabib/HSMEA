function srgtSRGT = srgtsRBNNFit(srgtOPT)
%Function srgtsRBNNFit fits the specified radial basis neural network model
%using the native MATLAB neural networks toolbox.
% 
%    srgtSRGT = srgtsRBNNFit(srgtOPT)
% 
%srgtSRGT is the surrogate structure that contains the following fields:
%* RBNN_NEWRBModel : radial basis neural network returned by the MATLAB
%                    function "newrb" (from the native MATLAB Neural
%                    Network toolbox).
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
%     options = srgtsRBNNSetOptions(X, Y);
% 
%     surrogate = srgtsRBNNFit(options)
% 
%     surrogate = 
% 
%     RBNN_NEWRBModel: [1x1 network]

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

srgtSRGT.RBNN_NEWRBModel = srgtsnewrb(...
    srgtOPT.P.', ...
    srgtOPT.T.', ...
    srgtOPT.RBNN_Goal, ...
    srgtOPT.RBNN_Spread, ...
    srgtOPT.RBNN_MN, ...
    srgtOPT.RBNN_DF, ...
    'off');

return
