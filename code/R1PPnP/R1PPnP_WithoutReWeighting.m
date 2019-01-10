%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the source code of R1PPnP without re-weighting
% We assume that this code runs outlier-free situation.
% <CopyRight 2017> Haoyin Zhou, Tao Zhang, Jayender Jagadeesan
% Please cite our paper when using this code:
% Haoyin Zhou, Tao Zhang, Jayender Jagadeesan, "Re-weighting and 1-Point RANSAC-Based PnP Solution to Handle Outliers", IEEE Transactions on Pattern Analysis and Machine Intelligence, 2018
%
% Pamameters:
% Xw: the three-dimensional world coordinates (3*N matrix, N is the number of points)
% uv2d: normalized image coordinate (2*N matrix, N is the number of points)
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ R_out, t_out] = R1PPnP_WithoutReWeighting( Xw, uv2d)

    N = size(Xw,2);
    uv = [uv2d;ones(1,N)];
 
    uv_centered = uv2d - mean(uv2d,2) * ones(1,N);
    uv_dis2center = sqrt(sum(uv_centered .* uv_centered));

    Xw_centered = Xw - mean(Xw,2) * ones(1,N);
    Xw_dis2center = sqrt(sum(Xw_centered .* Xw_centered));
    p_dis2center = uv_dis2center + 0.02 * Xw_dis2center;

    [uv_dis2center_min, min_idx] = min(p_dis2center);

    benchPid = min_idx;  
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

    uv_getridofbench_norm2vector = (sum(uv_getridofbench .* uv_getridofbench, 1));
            
    OneMatrix = ones(1, N-1);
    S = Xw_getridofbench - Xw(:,benchPid) * OneMatrix;

    uv_benchMatrix = [uv(:,benchPid)] * ones(1,N-1);
    
    uv_to_bench = uv_getridofbench(1:2,:) - uv_benchMatrix(1:2,:);
    
    R = eye(3); % initialiaztion
    mu = 0.001; % initialiaztion
    
    convertflag = 0;
    iter = 0;
    
    while(1)
        iter = iter + 1;
        p = mu * R * S + uv_benchMatrix;
        lambda = sum(uv_getridofbench .* p) ./ uv_getridofbench_norm2vector;

        if iter > 1
            if (norm(R-R_last) < 1e-4 && convertflag == 0 && det(R) == -1)
                lambda = 1.0 ./ lambda;
                convertflag = 1;
            end
        end
        lambdaMatrix = [lambda;lambda;lambda];
        q = lambdaMatrix .* uv_getridofbench;   
             
        S_q = (q - uv_benchMatrix);             
        
        lambdaMatrix_inv = 1.0 ./ lambdaMatrix;
        lambda_S = S .* lambdaMatrix_inv;
        lambda_S_q = S_q .* lambdaMatrix_inv;
        
        [UR,SR,VR] = svd(lambda_S_q * lambda_S');
        R = UR * VR'; 
        if (convertflag == 1 && det(R) == -1)
            R(3,:) = -R(3,:);
        end
        
        p_proj = [p(1,:)./p(3,:);p(2,:)./p(3,:)];
        P_proj_to_bench = p_proj - uv_benchMatrix(1:2,:);
        delta_mu = norm([uv_to_bench(1,:) uv_to_bench(2,:)]) / norm([P_proj_to_bench(1,:) P_proj_to_bench(2,:)]);
        mu = mu * delta_mu;

       if (iter == 1)
            done = false;
        else
            done = shouldbreak(iter, R, R_last);
        end
            
        if (done == true)
            break;
        end           

        R_last = R;
    end
    
       
   R_out = R;
   realq = 1 / mu * [1;1;1] * [1,lambda] .* [uv(:,benchPid),uv_getridofbench];
   t_raw = realq - R * [Xw(:,benchPid),Xw_getridofbench];
   t_out = sum(t_raw,2) / N;            

   [R_out, t_out] = finalrefinestep( R_out, t_out, Xw, uv, N);
       
end

function [done] = shouldbreak(iter, R, R_last)
   
    if (iter > 200)
        done = true;
        return;
    end    
  
    if (det(R) > 0.0 && norm(R - R_last) < 1e-4)
        done = true;
        return;
    end  
    
    done = false;
end









