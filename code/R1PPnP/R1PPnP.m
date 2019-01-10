%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full code of the R1PPnP algorithm
% <CopyRight 2017> Haoyin Zhou, Tao Zhang, Jayender Jagadeesan
%
% Please cite our paper when using this code:
% Haoyin Zhou, Tao Zhang, Jayender Jagadeesan, "Re-weighting and 1-Point RANSAC-Based PnP Solution to Handle Outliers", IEEE Transactions on Pattern Analysis and Machine Intelligence, 2018
%
% This algorithms randomly select one point as the benchmark and do
% iterative calculation.
% It is worth noting that this algorithm is not intended for small points sets,
% the number of inliiers should be larger than 20. 
%
% Input Pamameters:
% Xw: the three-dimensional world coordinates (3*N matrix, N is the number of points)
% uv2d: normalized image coordinate (2*N matrix, N is the number of points)
% focalLength: camera focalLength
% inlierErrorPixelTH: if the reporjection error is smaller than
%                     inlierErrorPixelTH, then it is considered as an
%                     inlier. 10 is the default value
% 
% Output Pamameters:
% R_out: rotation estimation result
% t_out: translation estimation result
% trialcount: Total number of RANSAC trials
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ R_out, t_out, trialcount] = R1PPnP( Xw, uv2d, focalLength, inlierErrorPixelTH)

    if (nargin < 4)
        inlierErrorPixelTH = 10;     
    end

    R_out = eye(3);
    t_out = [0;0;1];      
    
    N = size(Xw,2);
    uv = [uv2d;ones(1,N)];
    
    MaxTry = min(N,500);    % number of initial value try
 
    bestinliners_npts = 0;    
    
    inliners_npts_record = zeros(500,1);
  
    OneMatrix = ones(1, N-1);
 
    uv_centered = uv2d - mean(uv2d,2) * ones(1,N);
    uv_dis2center = sqrt(sum(uv_centered .* uv_centered));
    Xw_centered = Xw - mean(Xw,2) * ones(1,N);
    Xw_dis2center = sqrt(sum(Xw_centered .* Xw_centered));
    p_dis2center = uv_dis2center + 0.02 * Xw_dis2center;
    
     [sorted_uv_dis2center, sorted_idx] = sort(p_dis2center);
 
    for ransactry = 1:MaxTry
         benchPid = sorted_idx(ransactry);  % randomly select a point as the benchmark.
          
        if (benchPid == 1)
            Xw_getridofbench = [Xw(:,2:N)];
            uv_getridofbench = [uv(:,2:N)];
        elseif (benchPid == N)
            Xw_getridofbench = [Xw(:,1:N-1)];
            uv_getridofbench = [uv(:,1:N-1)];
        else
            Xw_getridofbench = [Xw(:,1:benchPid-1), Xw(:,benchPid+1:N)];
            uv_getridofbench = [uv(:,1:benchPid-1), uv(:,benchPid+1:N)];
        end
        uv_getridofbench_norm2vector = 1 ./ (sum(uv_getridofbench .* uv_getridofbench, 1));

        S = Xw_getridofbench - Xw(:,benchPid) * OneMatrix;

        uv_benchMatrix = [uv(:,benchPid)] * ones(1,N-1);
        uv_to_bench = uv_getridofbench(1:2,:) - uv_benchMatrix(1:2,:);

        R = eye(3); % initialiaztion
        mu = 0.001; % initialiaztion
      
        iter = 0;        
        while (1)
            iter = iter + 1;
            p = mu * R * S + uv_benchMatrix;
            lambda = sum(uv_getridofbench .* p) .* uv_getridofbench_norm2vector;
            lambdaMatrix = [lambda;lambda;lambda];

            q = lambdaMatrix .* uv_getridofbench;   

           lambdaMatrix_inv = 1.0 ./ lambdaMatrix;

            p_proj = [p(1,:)./p(3,:);p(2,:)./p(3,:)];   
            costvector = focalLength * (sum((p_proj - uv_getridofbench(1:2,:)).^2)).^0.5;

            weight = inlierErrorPixelTH ./ costvector;
            weight(weight > 1.0) = 1.0;
            weightMatrix = [1;1;1] * weight;  
            weightMatrix = weightMatrix .* lambdaMatrix_inv;
            
            wS_centered = weightMatrix .* S;
            wS_q = weightMatrix .* (q - uv_benchMatrix);              
 
            [UR,SR,VR] = svd(wS_q * wS_centered');
            R = UR * VR'; 

            P_proj_to_bench_weighted = (p_proj - uv_benchMatrix(1:2,:)) .* weightMatrix(1:2,:);
            uv_to_bench_weighted = uv_to_bench(1:2,:) .* weightMatrix(1:2,:);
            delta_mu = norm([uv_to_bench_weighted(1,:) uv_to_bench_weighted(2,:)]) / norm([P_proj_to_bench_weighted(1,:) P_proj_to_bench_weighted(2,:)]);
            mu = mu * delta_mu;
            
            inlinersidx = (costvector < inlierErrorPixelTH);
            inliners_npts = sum(inlinersidx);
            sumweight = sum(weight);
          
            if det(R) < 0
                R(3,:) = -R(3,:);
            end            
            
            if (iter > 29)
%                 fprintf('benchPid = %d, iter = %d, inliners_npts = %d (%d), sumweight = %.5f (%.5f), R_diff = %.5f\n', benchPid, iter,inliners_npts, inliners_npts - inliners_npts_record(iter - 20), sumweight, sumweight - sumweight_last, norm(R - R_last));
                done = shouldbreak(iter, norm(R-R_last), inliners_npts - inliners_npts_record(iter - 20));
                if (done == true)
                    break;
                end
            end
            
%             fprintf('iter = %d, inliners_npts = %d, sumweight = %.5f\n', iter, inliners_npts, sumweight);
            
            inliners_npts_record(iter) = inliners_npts;
            R_last = R;
            sumweight_last = sumweight;
         end
      
        %%
       inlinersidx = (costvector < inlierErrorPixelTH);
       inliners_npts = sum(inlinersidx);
       
       if (inliners_npts == 0)
           continue;
       end
       if (inliners_npts > bestinliners_npts)
           bestinliners_npts = inliners_npts;
           R_out = R;

           outliersidx = (inlinersidx == 0);
           uv_best = [uv(:,benchPid),uv_getridofbench(:,inlinersidx),uv_getridofbench(:,outliersidx)];
           Xw_best = [Xw(:,benchPid),Xw_getridofbench(:,inlinersidx),Xw_getridofbench(:,outliersidx)];
           realq = 1 / mu * [1;1;1] * [1,lambda(inlinersidx)] .* uv_best(:,1:(inliners_npts + 1));
           t_raw = realq - R * Xw_best(:,1:(inliners_npts + 1));
           t_out = sum(t_raw,2) / (inliners_npts + 1); 

           fracinliers = (inliners_npts+1) / N;
           fracinliers = max(eps, fracinliers);
           
           if (fracinliers > 0.55)
               punish = 0.0;
           else
               punish =  -4.6052 / log(1 - fracinliers);%-4.6052 / log(1 - fracinliers); 
           end
       
%            fprintf('try: benchPid = %d, iter = %d, inliners_ratio = %.1f, punish = %.1f\n', benchPid, iter, 100 * (inliners_npts+1) / N, punish);
       end
       if (punish < ransactry)
%             fprintf('break! punish = %f, ransactry = %d\n', punish, ransactry);
            [R_out, t_out] = finalrefinestep_weighted( R_out, t_out, Xw_best, uv_best, N, bestinliners_npts+1, focalLength, inlierErrorPixelTH);
            trialcount = ransactry;
            break;
       end
    end
    
end


 function [done] = shouldbreak(iter, R_diff, inliners_npts_diff)

    if (iter > 500)
        done = true;
        return;
    end
    
    if (inliners_npts_diff < 1 )
        done = true;
        return;
    end
 
    if (R_diff < 1e-4)
        done = true;
        return;
    end
    
    done = false;
 end








