function r = krig_corrgauss(theta,d)
% Gaussian correlation function
%===============================
%           n
%   r_i = prod exp(-theta_j * d_ij^2) ,  i = 1,...,m
%          j=1
%
% If length(theta) = 1, then the model is isotropic:
% all  theta_j = theta .
%
% Call:    r = dace_corrgauss(theta, d)
%          [r, dr] = dace_corrgauss(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation

[m, n] = size(d);  % number of differences and dimension of data
if  length(theta) == 1
    theta = repmat(theta,1,n);
elseif  length(theta) ~= n
    error(['Length of theta must be 1 or ',num2str(n)]);
end

td = d.^2 .* repmat(-theta(:).',m,1);
r = prod(exp(td), 2);
return
