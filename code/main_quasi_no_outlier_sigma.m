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
nls   = [0:0.5:10]; % level of noise
npts  = 100; % number of points
pouts = 0; % number of outliers
num   = 20; % total number of trials

% compared methods
A= zeros(size(npts));
B= zeros(num,1);

name= {'LHM', 'EPnP+GN', 'RPnP', 'DLS',          'OPnP', 'ASPnP', 'SDP',    'PPnP', 'EPPnP', 'REPPnP', 'R1PPnP'};
f= {    @LHM, @EPnP_GN,  @RPnP, @robust_dls_pnp, @OPnP, @ASPnP,   @GOP,     @PPnP,   @EPPnP,  @REPPnP, @R1PPnP_WithoutReWeighting};
f2 = {[],[],[],[],[],[],[],[],[],[],[]}; %post ransac method
ransacsamples = {0,0,0,0,0,0,0,0,0,0,0}; %number of samples to apply ransac
marker          = {'x',   's',      'd',      '^',    '>',   '<',      'v',      '*',     '+',         '^',    'x'};
color           = {'r',   'g',      [1,0.5,0],'m',    'c',   'b',      'y',      [1,0.5,1],[0,0.5,1],  'r',    'k'};
markerfacecolor=  {'r',   'g',      [1,0.5,0],'m',    'c',   'b',      'y',      'n',     'n',         'r',    'k'};


method_list= struct('name', name, 'f', f,'f2', f2, 'ransac',ransacsamples, 'mean_r', A, 'mean_t', A,...
    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, 'r', B, 't', B,...
    'marker', marker, 'color', color, 'markerfacecolor', markerfacecolor);

% experiments
for i= 1:length(nls)
       
        npt  = npts;
        pout = pouts;
        nl = nls(i);
        
        fprintf('npt = %d (sg = %d px) (%2.0f%%): ', npt, nl, pout*100);
    
        for k= 1:length(method_list)
            method_list(k).c = zeros(1,num);
            method_list(k).e = zeros(1,num);
            method_list(k).r = zeros(1,num);
            method_list(k).t = zeros(1,num);
        end

        %index_fail = [];
        index_fail = cell(1,length(name));
        
        for j= 1:num

            % camera's parameters
            width = 640;
            height = 480;
            f       = 1000;
            
            % generate 3d coordinates in camera space
            Xc  = [xrand(1,npt,[1 2]); xrand(1,npt,[1 2]); xrand(1,npt,[4 8])];
            t   = mean(Xc,2);
            R   = rodrigues(randn(3,1));
            XXw = inv(R)*(Xc-repmat(t,1,npt));

            % projection
            xx  = [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
            xxn = xx+randn(2,npt)*nl;
            xxn = xxn;

            %generate outliers (assigning to some 3D points more than one 2D correspondence)
            if (pout ~= 0)
            nout = max(1,round(npt * pout)); %at least we want 1 outlier
            idx  = randi(npt,1,nout);
            XXwo = XXw(:,idx);
            else
               nout = 0;
               XXwo = []; 
            end
            % assignation of random 2D correspondences
            xxo  = [xrand(1,nout,[-width/2 width/2]); xrand(1,nout,[-height/2 height/2])];
            
            labels = zeros(1,npt+nout);
            labels(1:npt) = 1;
            
            % pose estimation
            R1 = []; t1 = []; inliers = [];
            for k= 1:length(method_list)
                 if strcmp(method_list(k).name, 'R1PPnP')
                   tic;
                   [R1,t1]= method_list(k).f(XXw,xxn/f);
                   tcost = toc;
                 elseif strcmp(method_list(k).name, 'Reproj')
                    [R1,t1]= method_list(k).f([XXw, XXwo],[xxn, xxo]/f,R,t);
                 else
                    try
                        if (method_list(k).ransac == 0)
                            tic;
                            [R1,t1]= method_list(k).f([XXw, XXwo],[xxn, xxo]/f);
                            tcost = toc;
                        else
                            thr = 10; %ransac bandwidth in pixels
                            s = method_list(k).ransac;
                            X = [XXw, XXwo];
                            x = [xxn, xxo];
                            tic;
                            [R1, t1, inliers] = ransac(f, X,x, method_list(k).f, s, thr);
                            
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
            method_list(k).c(index_fail{k}) = [];
            method_list(k).e(index_fail{k}) = [];
            method_list(k).r(index_fail{k}) = [];
            method_list(k).t(index_fail{k}) = [];

            % computational cost should be computed in a separated procedure as
            % in main_time.m
            
            method_list(k).pfail(i) = 100 * numel(index_fail{k})/num;
            
            method_list(k).mean_c(i)= mean(method_list(k).c);
            method_list(k).mean_e(i)= mean(method_list(k).e);
            method_list(k).med_c(i)= median(method_list(k).c);
            method_list(k).med_e(i)= median(method_list(k).e);
            method_list(k).std_c(i)= std(method_list(k).c);
            method_list(k).std_e(i)= std(method_list(k).e);

            method_list(k).mean_r(i)= mean(method_list(k).r);
            method_list(k).mean_t(i)= mean(method_list(k).t);
            method_list(k).med_r(i)= median(method_list(k).r);
            method_list(k).med_t(i)= median(method_list(k).t);
            method_list(k).std_r(i)= std(method_list(k).r);
            method_list(k).std_t(i)= std(method_list(k).t);
           
            %results deleting solutions where not all the methods finds one
            tmethod_list.c(unique([index_fail{:}])) = [];
            tmethod_list.e(unique([index_fail{:}])) = [];
            tmethod_list.r(unique([index_fail{:}])) = [];
            tmethod_list.t(unique([index_fail{:}])) = [];
            
            method_list(k).deleted_mean_c(i)= mean(tmethod_list.c);
            method_list(k).deleted_mean_e(i)= mean(tmethod_list.e);
            method_list(k).deleted_med_c(i)= median(tmethod_list.c);
            method_list(k).deleted_med_e(i)= median(tmethod_list.e);
            method_list(k).deleted_std_c(i)= std(tmethod_list.c);
            method_list(k).deleted_std_e(i)= std(tmethod_list.e);

            method_list(k).deleted_mean_r(i)= mean(tmethod_list.r);
            method_list(k).deleted_mean_t(i)= mean(tmethod_list.t);
            method_list(k).deleted_med_r(i)= median(tmethod_list.r);
            method_list(k).deleted_med_t(i)= median(tmethod_list.t);
            method_list(k).deleted_std_r(i)= std(tmethod_list.r);
            method_list(k).deleted_std_t(i)= std(tmethod_list.t);
        end
       
end

save ../results/quasiresultsSigma method_list npt pout nls;

plotQuasisigmas;
saveas(1, '../results/figure9(a).png');
saveas(3, '../results/figure9(b).png');