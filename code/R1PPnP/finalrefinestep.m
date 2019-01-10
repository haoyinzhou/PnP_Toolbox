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

function [R, T] = finalrefinestep( R, T, Xw, uv, N)

    meanXw = mean(Xw,2);
    S = (Xw - meanXw * ones(1,N));
    uv_norm2vector = (sum(uv .* uv, 1));                
    
    for iter = 1:200
        p = R * Xw + T * ones(1,N);
        lambda = sum(uv .* p) ./ uv_norm2vector;
        lambda_matrix = [lambda;lambda;lambda];
        q = lambda_matrix .* uv;   
       
        meanq = mean(q,2);
        q_centered = q - meanq * ones(1,N);

        lambda_matrix_inv = 1.0 ./ lambda_matrix;
        lambda_S = S .* lambda_matrix_inv;
        lambda_q_centered = q_centered .* lambda_matrix_inv;
        
        [UR,SR,VR] = svd(lambda_q_centered * lambda_S');
        R = UR * VR'; 
        if (det(R) < 0)
            R(3,:) = -R(3,:);
        end  
        
        t_raw = q - R * Xw;        
        t_raw_weighted = t_raw .* lambda_matrix_inv;
        T = sum(t_raw_weighted,2) / sum(lambda_matrix_inv(1,:));
        
        if iter > 1
            Rdiff = norm(R - R_last);
            if Rdiff < 5e-5
                break;
            end
         end
         R_last = R;
     end
end





