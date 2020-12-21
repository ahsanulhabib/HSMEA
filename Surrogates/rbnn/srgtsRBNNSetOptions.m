function srgtOPT = srgtsRBNNSetOptions(P, T, RBNN_Goal, RBNN_Spread, RBNN_MN, RBNN_DF)
%Function srgtsRBNNSetOptions creates the SURROGATES Toolbox option
%structure for radial basis neural network models. This structure contains
%the following fields:
%
%* GENERAL PARAMETERS
%
%   SRGT - Identifier of the surrogate technique: 'RBNN'.
%   P    - NPOINTS-by-NDV matrix, where NPOINTS is the number of points of
%          the sample and NDV is the number of design variables.
%          Default: Empty matrix.
%   T    - NPOINTS-by-1 vector of responses on the P matrix points.
%          Default: Empty vector.
%
%* RADIAL BASIS NEURAL NETWORK PARAMETERS
%
%   RBNN_Goal   - Mean squared error goal [ scalar ].
%                 Default: (0.05*mean(T) )^2.
%   RBNN_Spread - Spread of radial basis functions [ positive scalar ].
%                 Default: 1.
%   RBNN_MN     - Maximum number of neurons [ positive integer ].
%                 Default: NPOINTS.
%   RBNN_DF     - Number of neurons to add between displays [ positive
%                 integer ]. Default: 1.
%
%The SURROGATES Toolbox uses the native MATLAB neural network toolbox.
%
%This is how you can use srgtsRBNNSetOptions:
% 
%     OPTIONS = srgtsRBNNSetOptions: creates a structure with the empty
%     parameters.
%
%     OPTIONS = srgtsRBNNSetOptions(P, T): Given the sampled data P (input
%     variables) and T (output variables), it creates a structure with
%     default parameters used for all not specified fields.
%
%     OPTIONS = srgtsRBNNSetOptions(P, T, ...
%     RBNN_Goal, RBNN_Spread, RBNN_MN, RBNN_DF): it creates a structure
%     with each of the specified fields.
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
%     options = srgtsRBNNSetOptions(X, Y)
%
%     options =
%
%            SRGT: 'RBNN'
%               P: [5x1 double]
%               T: [5x1 double]
%       RBNN_Goal: 1.0000e-004
%     RBNN_Spread: 1
%         RBNN_MN: 5
%         RBNN_DF: 1

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
srgtOPT.SRGT = 'RBNN';

srgtOPT.P = P;
srgtOPT.T = T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% svr options
if nargin == 2
    NbPoints = size(P,1);
    
    srgtOPT.RBNN_Goal   = (0.05*mean(T) )^2;
    srgtOPT.RBNN_Spread = 1;
    srgtOPT.RBNN_MN     = NbPoints;
    srgtOPT.RBNN_DF     = 1;
else
    srgtOPT.RBNN_Goal   = RBNN_Goal;
    srgtOPT.RBNN_Spread = RBNN_Spread;
    srgtOPT.RBNN_MN     = RBNN_MN;
    srgtOPT.RBNN_DF     = RBNN_DF;
end

return
