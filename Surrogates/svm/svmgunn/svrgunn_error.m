function err = svrgunn_error(trnX,tstX,tstY,ker,beta,bias,loss,e)
%SVRERROR Calculate SVR Error
%
%  Usage: err = svrgunn_error(trnX,tstX,tstY,ker,beta,bias,loss,e)
%
%  Parameters: trnX   - Training inputs
%              tstX   - Test inputs
%              tstY   - Test targets
%              ker    - kernel function
%              beta   - Difference of Lagrange Multipliers
%              bias   - bias
%              loss   - loss function
%              e      - e insensitive
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

n = size(trnX,1);
m = length(tstY);
H = zeros(m,n);
for i=1:m
    for j=1:n
        H(i,j) = svmgunn_kernelmatrix(ker,tstX(i,:),trnX(j,:));
    end
end

errvec = (H*beta - tstY + bias);

switch lower(loss)
    case 'einsensitive',
        errvec = abs(errvec) - e;
        err = mean(errvec.*(errvec > 0));
    case 'quadratic',
        err = mean(errvec.*errvec);
    otherwise, disp('Error: Unknown Loss Function\n');
end

return
