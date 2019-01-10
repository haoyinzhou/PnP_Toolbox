function plotPlanarModel(title,point31,point1,image1,v2d1,v3d,sf)

if ischar(image1)
    im1 = imread(image1);
    im1 = imresize(im1,sf);
else
    im1 = image1;
end

% Show a figure with lines joining the accepted matches.
figure('name',title);
colormap('gray');
subplot(1,2,1); hold on;
xlabel('x');ylabel('y');zlabel('z');

hold on; line([v3d(1,[1:4,1]) ],[v3d(2,[1:4,1]) ],[v3d(3,[1:4,1]) ],'color','b','LineWidth',4);

plot3(point31(1,:),point31(2,:),point31(3,:),'.g','MarkerSize',7);

axis square;
subplot(1,2,2);
imagesc(im1);
hold on;

%show vertex projection 
hold on; line([v2d1(1,1:4) v2d1(1,1)],[v2d1(2,1:4) v2d1(2,1)],'color','b','LineWidth',4);


hold on; plot(v2d1(1,1),v2d1(2,1),'o','MarkerSize',15,'MarkerFaceColor','c');
plot(v2d1(1,2),v2d1(2,2),'o','MarkerSize',15,'MarkerFaceColor','g');
plot(v2d1(1,4),v2d1(2,4),'o','MarkerSize',15,'MarkerFaceColor','y');

plot(point1(1,:),point1(2,:),'.g','MarkerSize',7);

axis off;
hold off;

end