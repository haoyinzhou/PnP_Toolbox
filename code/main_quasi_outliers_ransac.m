%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please cite our paper when using this code:
% Haoyin Zhou, Tao Zhang, Jayender Jagadeesan, "Re-weighting and 1-Point RANSAC-Based PnP Solution to Handle Outliers", IEEE Transactions on Pattern Analysis and Machine Intelligence, 2018
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
IniToolbox;

% experimental parameters
nl = 5; % level of noise
npts  = 100;  % number of points
pouts = [0.00:0.05:0.80];  % percentage of outliers
num   = 20; % total number of trials

% compared methods
A= zeros(size(npts));
B= zeros(num,1);

name= {'RNSC P3P','RNSC RP4P RPnP','RNSC P3P OPnP','RNSC P3P ASPnP', 'REPPnP', 'R1PPnP'};
f = {@kP3P,@RPnP,@kP3P,@kP3P,@REPPnP, @R1PPnP};
f2 = {[],@RPnP,@OPnP,@ASPnP,[],[]}; %post ransac method
ransacsamples = {3,4,3,3,0,0}; %number of samples to apply ransac
marker= {'o','d','>','<','^','x'};
color= {'g',[1,0.5,0],'c','b','r','k'};
markerfacecolor=  {'g',[1,0.5,0],'c','b','r','k'};

method_list= struct('name', name, 'f', f,'f2', f2, 'ransac',ransacsamples, 'mean_r', A, 'mean_t', A,...
    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, 'r', B, 't', B,...
    'marker', marker, 'color', color, 'markerfacecolor', markerfacecolor);

% experiments
for i= 1:length(pouts)
        
        npt  = npts;
        pout = pouts(i);
        
        fprintf('npt = %d (sg = %d px) (%2.0f%%): ', npt, nl, pout*100);
    
        for k= 1:length(method_list)
            method_list(k).c = zeros(1,num);
            method_list(k).e = zeros(1,num);
            method_list(k).r = zeros(1,num);
            method_list(k).t = zeros(1,num);
        end

        index_fail = cell(1,length(name));
        
        for j= 1:num

            % camera's parameters
            f  = 1000;
            
            % generate 3d coordinates in camera space
            Xc = [xrand(1,npt,[1 2]); xrand(1,npt,[1 2]); xrand(1,npt,[4 8])];
            t   = mean(Xc,2);
            R   = rodrigues(randn(3,1));
            XXw = inv(R)*(Xc-repmat(t,1,npt));
            
            % projection
            xx  = [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
           % randomvals = randn(2,npt);
            xxn = xx + randn(2,npt) * nl;
                         
            %generate outliers (assigning to some 3D points more than one 2D correspondence)
            if (pout ~= 0)
                nout = max(1,round((npt * pout)/(1-pout))); %at least we want 1 outlier

                idx  = randi(npt,1,nout);
                XXwo = XXw(:,idx);
            else
               nout = 0;
               XXwo = []; 
            end
            % assignation of random 2D correspondences
            xxo  = [xrand(1,nout,[min(xxn(1,:)) max(xxn(1,:))]); xrand(1,nout,[min(xxn(2,:)) max(xxn(2,:))])];
            
            % test
            % pose estimation
            R1 = []; t1 = []; inliers = [];
            for k= 1:length(method_list)
                 if strcmp(method_list(k).name, 'R1PPnP')
                    tic;
                    minInlierErrorinPixel = 10.0;
                    [R1,t1]= method_list(k).f([XXw, XXwo],[xxn, xxo]/f, f, minInlierErrorinPixel);
                    tcost = toc;
                 elseif strcmp(method_list(k).name, 'REPPnP')
                    tic;
                    [R1,t1,~,robustiters,terr]= method_list(k).f([XXw, XXwo],[xxn, xxo]/f);
                    tcost = toc;
                 else
                    try
                        if (method_list(k).ransac == 0)
                            tic;
                            robustiters = 0;
                            [R1,t1]= method_list(k).f([XXw, XXwo],[xxn, xxo]/f);
                            tcost = toc;
                        else
                            thr = 10;
                            s = method_list(k).ransac;
                            X = [XXw, XXwo];
                            x = [xxn, xxo];
                            K = [f 0 0;0 f 0; 0 0 1];
                            
                            tic;
                            [R1, t1, inliers, robustiters] = ransac(K, X, x, method_list(k).f, s, thr);
                            
                            if (~isempty(method_list(k).f2) && ~isempty(inliers))
                                [R1,t1]= method_list(k).f2(X(:,inliers),x(:,inliers)/f);
                            end
                            tcost = toc;
                            
                        end
                    catch
                    %   disp(['The solver - ',method_list(k).name,' - encounters internal errors!!!\n']);
                       index_fail{k} = [index_fail{k}, j];
                        continue;
                       % break;
                    end
                 end

                %no solution
                if size(t1,2) < 1
                   % disp(['The solver - ',method_list(k).name,' - returns no solution!!!\n']);
                    index_fail{k} = [index_fail{k}, j];
                    continue;
                   % break;
                elseif (sum(sum(sum(imag(R1).^2))>0) == size(R1,3) || sum(sum(imag(t1(:,:,1)).^2)>0) == size(t1,2))
                    index_fail{k} = [index_fail{k}, j];
                    continue;
                end
                %choose the solution with smallest error 
                error = inf;
                for jjj = 1:size(R1,3)
                    if (sum(sum(imag(R1(:,:,jjj)).^2)) + sum(imag(t1(:,jjj)).^2) > 0)
                        break
                    end            
                    
                    tempy = cal_pose_err([R1(:,:,jjj) t1(:,jjj)],[R t]);
                    if sum(tempy) < error
                        cost  = tcost;
                        %L2 error is computed without taing into account the outliers
                        ercorr= mean(sqrt(sum((R1(:,:,jjj) * XXw +  t1(:,jjj) * ones(1,npt) - Xc).^2)));
                        y     = tempy;
                        error = sum(tempy);
                    end
                end

                method_list(k).robustiters(j) = robustiters;
                method_list(k).c(j)= cost * 1000;
                method_list(k).e(j)= ercorr;
                method_list(k).r(j)= y(1);
                method_list(k).t(j)= y(2);
            end

            showpercent(j,num);
        end
        fprintf('\n');
    
        % save result
        for k= 1:length(method_list)
            
            %results without deleting solutions
            tmethod_list = method_list(k);
            method_list(k).robustiters(index_fail{k}) = [];
            method_list(k).c(index_fail{k}) = [];
            method_list(k).e(index_fail{k}) = [];
            method_list(k).r(index_fail{k}) = [];
            method_list(k).t(index_fail{k}) = [];

            method_list(k).pfail(i) = 100 * numel(index_fail{k})/num;
            
            method_list(k).mean_rit(i) = mean(method_list(k).robustiters);
            method_list(k).mean_c(i)= mean(method_list(k).c);
            method_list(k).mean_e(i)= mean(method_list(k).e);
            method_list(k).med_rit(i) = median(method_list(k).robustiters);
            method_list(k).med_c(i)= median(method_list(k).c);
            method_list(k).med_e(i)= median(method_list(k).e);
            method_list(k).std_rit(i) = std(method_list(k).robustiters);
            method_list(k).std_c(i)= std(method_list(k).c);
            method_list(k).std_e(i)= std(method_list(k).e);

            method_list(k).mean_r(i)= mean(method_list(k).r);
            method_list(k).mean_t(i)= mean(method_list(k).t);
            method_list(k).med_r(i)= median(method_list(k).r);
            method_list(k).med_t(i)= median(method_list(k).t);
            method_list(k).std_r(i)= std(method_list(k).r);
            method_list(k).std_t(i)= std(method_list(k).t);
           
            %results deleting solutions where not all the methods finds one
            tmethod_list.robustiters(unique([index_fail{:}])) = [];
            tmethod_list.c(unique([index_fail{:}])) = [];
            tmethod_list.e(unique([index_fail{:}])) = [];
            tmethod_list.r(unique([index_fail{:}])) = [];
            tmethod_list.t(unique([index_fail{:}])) = [];
            
            method_list(k).deleted_mean_rit(i) = mean(tmethod_list.robustiters);
            method_list(k).deleted_mean_c(i)= mean(tmethod_list.c);
            method_list(k).deleted_mean_e(i)= mean(tmethod_list.e);
            method_list(k).deleted_med_rit(i) = median(tmethod_list.robustiters);
            method_list(k).deleted_med_c(i)= median(tmethod_list.c);
            method_list(k).deleted_med_e(i)= median(tmethod_list.e);
            method_list(k).deleted_std_rit(i) = std(tmethod_list.robustiters);
            method_list(k).deleted_std_c(i)= std(tmethod_list.c);
            method_list(k).deleted_std_e(i)= std(tmethod_list.e);

            method_list(k).deleted_mean_r(i)= mean(tmethod_list.r);
            method_list(k).deleted_mean_t(i)= mean(tmethod_list.t);
            method_list(k).deleted_med_r(i)= median(tmethod_list.r);
            method_list(k).deleted_med_t(i)= median(tmethod_list.t);
            method_list(k).deleted_std_r(i)= std(tmethod_list.r);
            method_list(k).deleted_std_t(i)= std(tmethod_list.t);
        end
        
        prealout(i) = nout/npt; %exact percent of outliers    
        
end

save ../results/quasiresultsOutliersRansac method_list npt pouts prealout;

plotQuasioutliersRansac;

saveas(1, '../results/figure11(c).png');
saveas(3, '../results/figure11(d).png');
saveas(8, '../results/figure12(quasi_case_not_in_paper).png');
