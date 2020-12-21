function k = svmgunn_kernelmatrix(ker,kerOptions,u,v)
%SVKERNEL kernel for Support Vector Methods
%
%  Usage: k = prvt_svrgunn_kernelmatrix(ker,u,v)
%
%  Parameters: ker - kernel type
%              u,v - kernel arguments
%
%  Values for ker: 'Linear'                  -
%                  'Polynomial'              - p1 is degree of polynomial
%                  'GaussianRBF'             - p1 is width of rbfs (sigma)
%                  'MultiLayerPerceptron'    - p1 is scale, p2 is offset
%                  'LinearSpline'            -
%                  'LinearBSpline'           - p1 is degree of bspline
%                  'TrigonometricPolynomial' - p1 is degree
%                  'ExponentialRBF'          - p1 is width of rbfs (sigma)
%                  'Anova(xxx)'              - p1 is max order of terms
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

switch ker
    case 'Linear'
        k = u*v';

    case 'Polynomial'
        p1 = kerOptions;
        k = (u*v' + 1)^p1;

    case 'GaussianRBF'
        p1 = kerOptions;
        k = exp(-(u-v)*(u-v)'/(2*p1^2));

    case 'ExponentialRBF'
        p1 = kerOptions;
        k = exp(-sqrt((u-v)*(u-v)')/(2*p1^2));

    case 'MultiLayerPerceptron'
        p1 = kerOptions(1);
        p2 = kerOptions(2);
        k = tanh(p1*u*v'/length(u) + p2);

    case 'TrigonometricPolynomial'
        p1 = kerOptions;
        z = sin(p1 + 1/2)*2*ones(length(u),1);
        i = find(u-v);
        z(i) = sin(p1 + 1/2)*(u(i)-v(i))./sin((u(i)-v(i))/2);
        k = prod(z);

    case 'LinearSpline'
        z = 1 + u.*v + (1/2)*u.*v.*min(u,v) - (1/6)*(min(u,v)).^3;
        k = prod(z);

    case 'LinearBSpline'
        p1 = kerOptions;
        z = 0;
        for r = 0: 2*(p1+1)
            z = z + (-1)^r*binomial(2*(p1+1),r)*(max(0,u-v + p1+1 - r)).^(2*p1 + 1);
        end
        k = prod(z);

    case 'AnovaSpline-1'
        z = 1 + u.*v + u.*v.*min(u,v) - ((u+v)/2).*(min(u,v)).^2 + (1/3)*(min(u,v)).^3;
        k = prod(z);

    case 'AnovaSpline-2'
        z = 1 + u.*v + (u.*v).^2 + (u.*v).^2.*min(u,v) - u.*v.*(u+v).*(min(u,v)).^2 + (1/3)*(u.^2 + 4*u.*v + v.^2).*(min(u,v)).^3 - (1/2)*(u+v).*(min(u,v)).^4 + (1/5)*(min(u,v)).^5;
        k = prod(z);

    case 'AnovaSpline-3'
        z = 1 + u.*v + (u.*v).^2 + (u.*v).^3 + (u.*v).^3.*min(u,v) - (3/2)*(u.*v).^2.*(u+v).*(min(u,v)).^2 + u.*v.*(u.^2 + 3*u.*v + v.^2).*(min(u,v)).^3 - (1/4)*(u.^3 + 9*u.^2.*v + 9*u.*v.^2 + v.^3).*(min(u,v)).^4 + (3/5)*(u.^2 + 3*u.*v + v.^2).*(min(u,v)).^5 - (1/2)*(u+v).*(min(u,v)).^6 + (1/7)*(min(u,v)).^7;
        k = prod(z);

    case 'AnovaBSpline'
        p1 = kerOptions;
        z = 0;
        for r = 0: 2*(p1+1)
            z = z + (-1)^r*binomial(2*(p1+1),r)*(max(0,u-v + p1+1 - r)).^(2*p1 + 1);
        end
        k = prod(1 + z);

    otherwise
        k = u*v';
end

return

function b = binomial(n,k)
%BINOMIAL compute binomial coefficient
%
%  Usage: b = binomial(n,k)
%
%  Parameters:     ( n )
%              b = (   )
%                  ( k )
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

  if (nargin ~=2 | k>n | k<0 | n<0) % check correct number of arguments
    help binomial
  else
     
     b = prod(k+1:n)/prod(1:(n-k));
     
  end
