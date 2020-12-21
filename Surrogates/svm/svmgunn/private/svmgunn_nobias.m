function nb = svmgunn_nobias(SVR_Kernel)
%NOBIAS returns true if SVM kernel has no implicit SVR_Bias
%
%  Usage: nb = svrgunn_nobias(SVR_Kernel)
%
%  Parameters: SVR_Kernel - kernel type
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

switch lower(SVR_Kernel)
    case {'linear','sigmoid'}
        %,'anovaspline1','anovaspline2','anovaspline3'}
        nb = 1;
    otherwise
        nb = 0;
end

return
