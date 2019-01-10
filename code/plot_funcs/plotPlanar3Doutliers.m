load planar3DresultsOutliers

close all;
yrange= [0 2];

i= 0; w= 300; h= 300;

yrange= [0 max([method_list(:).deleted_mean_r])];

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_mean_r','Mean Rotation Error',...
    '% of outliers','Rotation Error (degrees)');

figure('color','w','position',[w*i,100+h,w,h]);i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_med_r','Median Rotation Error',...
    '% of outliers','Rotation Error (degrees)');

yrange= [0 max([method_list(:).deleted_mean_t])];

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_mean_t','Mean Translation Error',...
    '% of outliers','Translation Error (%)');

figure('color','w','position',[w*i,100+h,w,h]);i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_med_t','Median Translation Error',...
    '% of outliers','Translation Error (%)');

yrange= [0 max([method_list(:).deleted_mean_e])];

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_mean_e','Mean L2 Error',...
    '% of outliers','L2 error');

figure('color','w','position',[w*i,100+h,w,h]);
xdrawgraph(pouts*100,yrange,method_list,'deleted_med_e','Median L2 Error',...
    '% of outliers','L2 error');

yrange= [0 min(max(1,2*max([method_list(:).pfail])),100)];
figure('color','w','position',[w*i,100+2*h,w,h]);%i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'pfail','No solution x method',...
    '% of outliers','% method fails');
i=i+1;

% yrange= [0 2*max([method_list(:).med_e])];
% 
% figure('color','w','position',[w*i,100,w,h]);%i=i+1;
% xdrawgraph(pouts*100,yrange,method_list,'deleted_mean_e','Mean L2 Error',...
%     '% of outliers/inliers','L2 error');
% 
% figure('color','w','position',[w*i,100+h,w,h]);i=i+1;
% xdrawgraph(pouts*100,yrange,method_list,'deleted_med_e','Median L2 Error',...
%     '% of outliers/inliers','L2 error');
%  

yrange= [0 max([method_list(:).deleted_mean_c])];

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_mean_c','Mean Cost',...
    '% of outliers','Cost (ms)');

figure('color','w','position',[w*i,100+h,w,h]);i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_med_c','Median Cost',...
    '% of outliers','Cost (ms)');
