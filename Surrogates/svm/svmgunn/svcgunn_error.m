function err = svcerror(trnX,trnY,tstX,tstY,ker,alpha,bias)
%SVCERROR Calculate SVC Error
%
%  Usage: err = svcerror(trnX,trnY,tstX,tstY,ker,alpha,bias)
%
%  Parameters: trnX   - Training inputs
%              trnY   - Training targets
%              tstX   - Test inputs
%              tstY   - Test targets
%              ker    - kernel function
%              beta   - Lagrange Multipliers
%              bias   - bias
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

  if (nargin ~= 7) % check correct number of arguments
    help svcerror
  else

    n = size(trnX,1);
    m = length(tstY);
    H = zeros(m,n);  
    for i=1:m
      for j=1:n
        H(i,j) = trnY(j)*svkernel(ker,tstX(i,:),trnX(j,:));
      end
    end
    predictedY = sign(H*alpha + bias);
    err = sum(predictedY ~= tstY);

  end
