function r = krig_matern_5_2(theta,d)
% Matern 5/2 correlation function
%=================================
%          d
%   r_i = prod (1 + sqrt(5)*theta_l*m_i + 5/3 theta_l^2*m_i^2 ) * exp(-sqrt(5)*theta_l*m_il) ;  
%         j=1                                                             where, i = 1,...,n, 
%                                                                            and j = 1,...,d
%          
%
%   m_i = |(x_i_j - x'_i_j)|
%
% If length(theta) = 1, then the model is isotropic, i.e. all  theta_j = theta
%
% Input:
% -------------------------------------------------------------------------------------------------
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points weighted
%          with Wstar obtained from ECR for each Principal Components 
% Output:
% -------------------------------------------------------------------------------------------------
% r     :  correlation value for a particular Principal Component
n = size(d,1);
Theta = theta(:).';
term1 = (1 + sqrt(5)*repmat(Theta,n,1).*abs(d) + 5/3*(repmat(Theta,n,1)).^2.*(d.^2)); 
term2 = exp(-sqrt(5)*repmat(Theta,n,1).*abs(d));
r = prod(term1.*term2,2);
return
