function [sumArea2,sumArea]=labeledImageArea(imMedian,labeledImage,centroids,foldername,iprob,N)
%%%%%%%%%%%
%
%%%%%%%%%%
% coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels
figure;imagesc(imMedian);colormap(gray);hold on; axis off;contour(labeledImage,'yo')
%s  = regionprops(labeledImage, 'centroid');
%figure;imagesc(labeledImage);axis off
hold on
s  = regionprops(labeledImage, 'Area');
Area = cat(1, s.Area);
text(centroids(:,1)-3, centroids(:,2), num2str(Area),'FontSize',10,'Color','r')
hold off
%figure; imagesc(coloredLabels);colormap(gray), axis off
s=sprintf('print -depsc %s/Area2_%d_Ny%d,print -djpeg %s/Area2_%d_Ny%d;',foldername,iprob,N,foldername,iprob,N); eval(s)
 
% figure;imagesc(labeledImage);axis off
% hold on
figure;imagesc(imMedian);colormap(gray);hold on; axis off;contour(labeledImage,'yo')
text(centroids(:,1)-3, centroids(:,2), num2str(round(Area*2.*(0.138)^2)),'FontSize',10,'Color','r')
hold off
s=sprintf('print -depsc %s/Area_%d_Ny%d,print -djpeg %s/Area_%d_Ny%d;',foldername,iprob,N,foldername,iprob,N); eval(s)
  

% figure;imagesc(coloredLabels);hold on; axis off;text(centroids(:,1)-3, centroids(:,2), num2str(round(Area*2.*(0.138)^2)),'FontSize',8,'Color',[1,0,1],'LineWidth',2)
% hold off
% s=sprintf('print -depsc %s/Area_colLab_%d_Ny%d,print -djpeg %s/Area_colLab_%d_Ny%d;',foldername,iprob,N,foldername,iprob,N); eval(s)
  sumArea2=sum(sum(Area(:)));
  AA2=Area*2.*(0.138)^2;
  sumArea=sum(sum(AA2(:)));


