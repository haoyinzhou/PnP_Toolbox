%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the final refine step of R1PPnP
% this step will slightly refine the results.
% <CopyRight 2017> Haoyin Zhou, Tao Zhang, Jayender Jagadeesan
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, T] = finalrefinestep_weighted( R, T, Xw, uv, N, Ninner0, focalLength, inlierErrorPixelTH)

    weight = [ones(1,Ninner0), zeros(1, N-Ninner0)];

    uv_norm2vector = (sum(uv .* uv, 1));
    meanXwInner = mean(Xw(:,1:Ninner0),2);
    S = Xw - meanXwInner * ones(1,N);
  
    for iter = 1:500
        p = R * Xw + T * ones(1,N);
        lambda = sum(uv .* p) ./ uv_norm2vector;
        lambda_matrix = [lambda;lambda;lambda];
        q = lambda_matrix .* uv;   
 
        lambda_matrix_inv = 1.0 ./ lambda_matrix;

        p_proj = [p(1,:)./p(3,:);p(2,:)./p(3,:)];   
        costvector = focalLength * (sum((p_proj - uv(1:2,:)).^2)).^0.5;
        weight(Ninner0+1:end) = (costvector(Ninner0+1:end) < inlierErrorPixelTH);
        weightMatrix = [1;1;1] * weight;                
       
        weightMatrix = weightMatrix .* lambda_matrix_inv;
        
        wS = weightMatrix .* S;
        
        meanq = mean(q(:,1:Ninner0),2);
        wSq = weightMatrix .* (q - meanq * ones(1,N));

        [UR,SR,VR] = svd(wSq * wS');
        R = UR * VR'; 
        if (det(R) < 0)
            R(3,:) = -R(3,:);
        end  
        
        t_raw = weightMatrix .* (q - R * Xw);
        T = sum(t_raw,2) ./ sum(weightMatrix,2);

        if iter > 1
            Rdiff = norm(R - R_last);
            if Rdiff < 5e-5
                break;
            end
         end
         R_last = R;
    end
end





