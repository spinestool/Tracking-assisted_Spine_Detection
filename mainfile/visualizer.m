
function retVal = visualizer(im1, im2, alphaVal, method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualizer for spine detection and segmentation
%
% author : M. Yagci & E.Erdil
% date   : 07.04.2012
%

if strcmp(method, 'RGB')    
    figure, imshow(im1), hold on
    hImage = imshow(im2); % TEST
    set(hImage, 'AlphaData', alphaVal);
else
    % non-rgb for uint8
    imshow(im1, 'DisplayRange', []), hold on
    colorizedImage = repmat(double(im2),[1 1 3]);
    colorizedImage(:,:,2,1) = 0 ; colorizedImage(:,:,3,1) = 0;
    hImage = imshow(colorizedImage);
    set(hImage, 'AlphaData', alphaVal);
end

