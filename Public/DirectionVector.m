function [w_all] = DirectionVector(M,pall)
count = 1;
w_all = [];
for k = 1:numel(pall)
    p = pall(k);
    NumPoints=zeros(1,M-1);
    for i=1:(M-1)
        NumPoints(i) = p;
    end
    % Partition the first objective
    Beta1 = [0:(NumPoints(1))]'/(NumPoints(1));
    % Save the previous values
    Beta = Beta1;
    for i = 2:M-1
        % Compute the combination i.e. p*(p-1)*(p-2)*,-----,*0
        ki = round((1-sum(Beta1,2))*(NumPoints(i)));
        Beta = [];
        for j =1:size(Beta1,1)
            % Compute each subvector of (0,1,...p)/p,(0,1,...p-1)/p,...
            BetaVec = [0:ki(j)]'/(NumPoints(i));
            numOfreplications = length(BetaVec); % identify the length
            % Replicate each of the previous values in the equal size to all the subvectors
            Beta = [Beta; [repmat(Beta1(j,:), numOfreplications,1) BetaVec] ];
        end
        Beta1 = Beta;
    end
    % Compute the last objective values
    BetaVec = 1 - sum(Beta1,2);
    w = [Beta BetaVec];%include the last objective values
    w = w*count + (1-count)/M;
    w_all = [w_all;w];
    count = count/2;
end
w_all(isnan(w_all(:,1)),:) = [];
end