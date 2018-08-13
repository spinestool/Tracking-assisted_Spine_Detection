function [waterShedSegmentation_all,L,waterShedLevels]=watershed_seg(Img,maxNumberOfPixels)

imEq = adapthisteq(Img); %Contrast-limited adaptive histogram equalization
%region based method that connects region with same intensity
imEMaxima = imextendedmax(Img,5);%computes the extended-maxima transform, which is the regional maxima of the H-maxima transform
imEMaxima = imfill(imEMaxima, 'holes');%fills holes in the binary image BW. 

imBW = im2bw(imEq, graythresh(imEq));%Convert an image to a binary image, based on threshold

% roi selection 
imBW = imdilate(imBW | imEMaxima, ones(3,3));

% find boundaries
imBW = imfill(imBW,'holes');

% watersheds (segmentation phase 1)
imEqComp = imcomplement(imEq);%%%Complement image.
% make background and imEMaxima ==> only maxima; modifies the intensity imEqComp using morphological reconstruction so it only has regional minima wherever imBW is nonzero.
imEqCompModified = imimposemin(imEqComp, ~imBW | imEMaxima);%%%%%Impose minima.
% apply watershed transform.
L = watershed(imEqCompModified); %figure,imagesc(L)

% graph-based (segmentation phase 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberOfComponents = max(max(L));
% maxNumberOfPixels = 6000;%%%800 better then 8000 on some cases ex 26. 
minNumberOfPixels = 20;  
nonzeroL = [];
waterShedSegmentation_all=zeros(size(Img));
for i = 1:numberOfComponents % i=0 is the watershed, ignore
%nnz - number of nonzero elements
    if((nnz(L == i) < maxNumberOfPixels) && (nnz(L == i) > minNumberOfPixels))
            imTempWaterShed1 = zeros(size(Img));
            waterShedSegmentation = (L == i);
            imTempWaterShed1(find(waterShedSegmentation == 1)) = 255;
            waterShedLevels(i).h = imTempWaterShed1; 
            waterShedSegmentation_all=waterShedSegmentation_all+waterShedSegmentation;

    end
end
