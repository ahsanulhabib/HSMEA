function dmodel = dace_build(X, Y)
% DACEMODEL - To construct kriging model
%
% Call
%    rmodel = dace_model(X, Y)
%
% Input
% X  : Data Points X(i,:), i=1,...,m
% 
% Output
% dmodel : DACE Model
%

% Check arguments
if nargin ~= 2
	error('dace_model requires 2 input arguments')
end

% Check design points
[m1, nx] = size(X);
[m2, ~] = size(Y);
if m1 ~= m2
	error('X and Y must have the same number of rows')
end

% theta = 10 * ones(1, nx);
% range = minmax(X');

dmodel = dacefit(X, Y, @regpoly1, @corrgauss, 10*ones(1,nx), ones(1,nx), 20*ones(1,nx));
% [dmodel, perf] = dacefit(X, Y, @regpoly1, @corrgauss, theta);
return
