function [fout] = NormDenorm(f,option,range)
if option == 1 % 1 = normalize;
%     try
%         if nargin == 3
%             fout = (f - repmat(range(1,:),size(f,1),1))./repmat((range(2,:)-range(1,:)),size(f,1),1);
%         else
%             fout = (f - repmat(min(f,[],1),size(f,1),1))./repmat((max(f,[],1)-min(f,[],1)),size(f,1),1);
%         end
%     catch
        if nargin == 3
            fout = zeros(size(f));
            minpoint = range(1,:);
            maxpoint = range(2,:);
            idscale = minpoint < maxpoint;
            N1 = size(f,1);
            if sum(idscale) > 0
                fout(:,idscale)=(f(:,idscale)-ones(N1,1)*minpoint(:,idscale))./(ones(N1,1)*maxpoint(:,idscale)-ones(N1,1)*minpoint(:,idscale));
            end
        else
            fout = zeros(size(f));
            minpoint = min(f,[],1);
            maxpoint = max(f,[],1);
            idscale = minpoint < maxpoint;
            N1 = size(f,1);
            if sum(idscale) > 0
                fout(:,idscale)=(f(:,idscale)-ones(N1,1)*minpoint(:,idscale))./(ones(N1,1)*maxpoint(:,idscale)-ones(N1,1)*minpoint(:,idscale));
            end
        end
%     end
else % 0 = denormalize
%     try
%         if nargin == 3
%             fout = repmat(range(1,:),size(f,1),1) + f.*repmat((range(2,:)-range(1,:)),size(f,1),1);
%         else
%             fout = repmat(min(f),size(f,1),1) + f.*repmat((max(f)-min(f)),size(f,1),1);
%         end
%     catch
        if nargin == 3
            fout = zeros(size(f));
            minpoint = range(1,:);
            maxpoint = range(2,:);
            idscale = minpoint < maxpoint;
            N1 = size(f,1);
            if sum(idscale) > 0
                fout(:,idscale)=repmat(minpoint(:,idscale),size(f,1),1) + f.*repmat((maxpoint(:,idscale)-minpoint(:,idscale)),size(f,1),1);
            end
        else
            fout = zeros(size(f));
            minpoint = min(f,[],1);
            maxpoint = max(f,[],1);
            idscale = minpoint < maxpoint;
            N1 = size(f,1);
            if sum(idscale) > 0
                fout(:,idscale)=repmat(minpoint(:,idscale),size(f,1),1) + f.*repmat((maxpoint(:,idscale)-minpoint(:,idscale)),size(f,1),1);
            end
        end
%     end
end
end