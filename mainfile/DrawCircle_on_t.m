
function [ph1,ph10]= DrawCircle_on_t( zd,terminating_pts1,r,r2)

%rad(1 : size(terminating_pts1(:, 1)))=5;
% for k = 1 : size(terminating_pts1, 1),
%      DrawCircle(terminating_pts1(k,1), terminating_pts1(k,2), rad(k), 32, 'b-');
%  end

[n,m]=size(zd);

 [ph1,ph10]=  DrawCircle1(n,m,terminating_pts1(:,2), terminating_pts1(:,1), r,r2);
 

% figure,imagesc(ph1),colormap(gray)
% hold on
% contour(ph1,[0 0], 'yo','LineWidth',2);  colormap(gray);hold off;

% figure,imagesc(ph10),colormap(gray)
% hold on
% contour(ph10,[0 0], 'yo','LineWidth',2);  colormap(gray);hold off;