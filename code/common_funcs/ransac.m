% This function is an optimization of the ransac function provided by:
% Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk

function [R, T, inliers, trialcount] = ransac(K, X, x, fittingfn, s, t, idsampl, maxTrials)

    % Test number of parameters
    narginchk ( 5, 8 ) ;
    
    [~, npts] = size(x);
     
    if nargin < 7; idsampl = 1:npts;    end;
    if nargin < 8; maxTrials = 5000;    end;
    
    p = 0.99;         % Desired probability of choosing at least one sample
                      % free from outliers

    bestR = NaN;      % Sentinel value allowing detection of solution failure.
    bestT = NaN;
    trialcount = 0;
    bestscore =  0;
    
    %we assume 50% of outliers (Dummy initialisation for number of trials.)
    N = log(1-p)/log(1-0.5^s); 
    nsampl = length(idsampl);
    
    while N > trialcount
      
        % Generate s random indicies in the range 1..npts
        ind = randsample(nsampl, s);
        ind = idsampl(ind);
        
        try
            %K = [f 0 0;0 f 0; 0 0 1];
            u = K \ [x(:,ind);ones(1,numel(ind))];
            u = u(1:2,:);
            [R, T] = fittingfn(X(:,ind), u);
        catch
            R = [];
            T = [];
        end

        [inliers, R, T] = evalRT(K, R, T, X, x, t);

        ninliers = length(inliers);
       
        if ninliers > bestscore    % Largest set of inliers so far...
            bestscore = ninliers;  % Record data for this model
            bestinliers = inliers;
            bestR = R;
            bestT = T;
            
            % Update estimate of N, the number of trials to ensure we pick,
            % with probability p, a data set with no outliers.
            fracinliers =  ninliers/nsampl;
            pNoOutliers = 1 -  fracinliers^s;
            pNoOutliers = max(eps, pNoOutliers);  % Avoid division by -Inf
            pNoOutliers = min(1-eps, pNoOutliers);% Avoid division by 0.
            N = log(1-p)/log(pNoOutliers);
            
%             fprintf('update N = %.3f, pNoOutliers = %.3f\n', N, pNoOutliers);
            
        end
        
%         fprintf('trialcount = %d, ninliers = %d, N = %.3f\n', trialcount, ninliers, N);
       
        trialcount = trialcount+1;
        
        % Safeguard against being stuck in this loop forever
        if trialcount > maxTrials
            break
        end
    end
    
    if ~isnan(bestR)   % We got a solution
        R = bestR;
        T = bestT;
        inliers = bestinliers;
    else
        R = [];
        T = [];
        inliers = [];
    end
end

function  [inliers, R, T] = evalRT(K, tR, tT, X, x, t)

R = [];
T = [];
inliers = [];

try
%K = [f 0 0;0 f 0; 0 0 1];
for i = 1:size(tR,3)
    newX  = K * (tR(:,:,i) * X + tT(:,i) * ones(1,size(X,2)));
    error = sqrt(sum((x - newX(1:2,:)./repmat(newX(3,:),2,1)).^2));
    tinliers = find(error < t);
    if length(tinliers) > length(inliers)
        inliers = tinliers;
        R = tR(:,:,i);
        T = tT(:,i);
    end
end
catch
    R = [];
    T = [];
    inliers = [];
end

end
