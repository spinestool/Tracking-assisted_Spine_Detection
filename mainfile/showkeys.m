% showkeys(image, locs)
%
% This function displays an image with SIFT keypoints overlayed.
%   Input parameters:
%     image: the file name for the image (grayscale)
%     locs: matrix in which each row gives a keypoint location (row,
%           column, scale, orientation)

function [filteredImage,filBinImg1]=showkeys(image, locs, group,iprob,subfoldername,N)
global rad rad2
disp('Drawing SIFT keypoints ...');

% Draw image with keypoints
figure,%('Position', [0 0 size(image,2) size(image,1)]);
colormap('gray');axis off
imagesc(image);
title(['SIFT keypoints for time series data' num2str(iprob)' ]);hold on;locs_spine=[];
% imsize = size(image);
for i = 1: size(locs,1)
    % Draw an arrow, each line transformed according to keypoint parameters.
    if(group(i) == 1)
        plot(locs(i, 2), locs(i, 1), 'rx', 'MarkerSize', 10);
        locs_spine=[locs_spine;locs(i, :)];
    else
        plot(locs(i, 2), locs(i, 1), 'bo', 'MarkerSize', 10);
    end
%     TransformLine(imsize, locs(i,:), 0.0, 0.0, 1.0, 0.0);
%     TransformLine(imsize, locs(i,:), 0.85, 0.1, 1.0, 0.0);
%     TransformLine(imsize, locs(i,:), 0.85, -0.1, 1.0, 0.0);
end
axis off
hold off;
 s=sprintf('print -depsc %s/nsift_cl_%d_N%d,print -djpeg %s/nsift_cl_%d_N%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)
%%%%%%%%%%%%%%%%added for the domain

[ph1,ph10]=DrawCircle_on_t( image,locs_spine,rad,rad2);
% keyboard
% zd_new=image.*ph1;
% figure, imagesc(image),hold on; colormap(gray);contour(zd_new,[0 0], 'ro');
% s=sprintf('print -dpng %s/cont_%d_N%d,print -depsc %s/cont_%d_N%d,print -djpeg %s/cont_%d_N%d;',subfoldername,iprob,N,subfoldername,iprob,N,subfoldername,iprob,N); eval(s)
% 
% filteredImage = bwmorph(zd_new, 'open');

filteredImage=bwmorph(ph1, 'open');
%figure, imagesc(filteredImage);colormap(gray)
% s=sprintf('print -depsc %s/bwmorph_%d_N%d,print -djpeg %s/bwmorph_%d_N%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)

filBinImg1=bwmorph(ph10, 'open');
%figure, imagesc(filBinImg1);colormap(gray)
% s=sprintf('print -depsc %s/bwmorph10_%d_N%d,print -djpeg %s/bwmorph10_%d_N%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)

%%%%%%%%%%%%%%%%

% ------ Subroutine: TransformLine -------
% Draw the given line in the image, but first translate, rotate, and
% scale according to the keypoint parameters.
%
% Parameters:
%   Arrays:
%    imsize = [rows columns] of image
%    keypoint = [subpixel_row subpixel_column scale orientation]
%
%   Scalars:
%    x1, y1; begining of vector
%    x2, y2; ending of vector
function TransformLine(imsize, keypoint, x1, y1, x2, y2)

% The scaling of the unit length arrow is set to approximately the radius
%   of the region used to compute the keypoint descriptor.
len = 6 * keypoint(3);

% Rotate the keypoints by 'ori' = keypoint(4)
s = sin(keypoint(4));
c = cos(keypoint(4));

% Apply transform
r1 = keypoint(1) - len * (c * y1 + s * x1);
c1 = keypoint(2) + len * (- s * y1 + c * x1);
r2 = keypoint(1) - len * (c * y2 + s * x2);
c2 = keypoint(2) + len * (- s * y2 + c * x2);

line([c1 c2], [r1 r2], 'Color', 'c');

