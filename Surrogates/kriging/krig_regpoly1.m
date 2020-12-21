function  [f, df] = krig_regpoly1(S)
%First order polynomial regression function
%
% Call:    f = dace_regpoly1(S)
%          [f, df] = dace_regpoly1(S)
%
% S : m*n matrix with design sites
% f = [1  s]
% df : Jacobian at the first point (first row in S) 

% hbn@imm.dtu.dk  
% Last update April 12, 2002

[m n] = size(S);
f = [ones(m,1)  S];

if  nargout > 1
  df = [zeros(n,1) eye(n)];
end

return
