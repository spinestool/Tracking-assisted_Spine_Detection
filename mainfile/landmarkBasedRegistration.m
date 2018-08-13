
function [bestTransformation,imreg] = landmarkBasedRegistration(im1, im2,desc1,locs1,desc2,locs2,foldernameout,step,iprob)
global registration_plot 

% [~, desc1, locs1] = sift(im1);
% [~, desc2, locs2] = sift(im2);
%%%registration visualization
if registration_plot==1
 figure; visualizer(uint8(im2), im2bw(im1, graythresh(im1)), 0.4, '');title('Before Registration')
s=sprintf('print -depsc %s/Before_reg%d_im_no%d,print -djpeg %s/Before_reg%d_im_no%d;',foldernameout,iprob,step,foldernameout,iprob,step); eval(s)
end
[match1, angles, locations] = match(im1, im2, desc1, desc2, locs1, locs2, 0.8,foldernameout,iprob,step);
indices = find(match1 ~= 0); %number of common key locs (like common spines but not nesesary)
match1 = match1(indices);%takes only matched points
angles = angles(indices); %takes angles of these matched lines (orientation)

roundedAngles = round(angles);
maxRepeatedAngles = mode(roundedAngles); %most repeated angles

indices = find(roundedAngles == maxRepeatedAngles);

translations = zeros(length(indices), 2); %(x, y)
for i = 1:length(indices)
    translations(i, 1) = locations(indices(i), 2) - size(im1, 2) - locations(indices(i), 1);
    translations(i, 2) = locations(indices(i), 4) - locations(indices(i), 3);
end

x = sort(translations(:, 1));
y = sort(translations(:, 2));
len = round(length(x) * 0.10);

x = x(len:(end - len)); 
y = y(len:(end - len));%%% we get as x and y ony 80 % of the elements
 
x = round(mean(x));
y = round(mean(y));

im1 = uint8(im1);
% visualizer(uint8(im2), im2bw(im1, graythresh(im1)), 0.4, '');
% s=sprintf('print -depsc %s/nonreg%d_im_no%d,print -djpeg %s/nonreg%d_im_no%d;',foldernameout,iprob,step,foldernameout,iprob,step); eval(s)

%translates im1 down/up by y axis and right/left by x axis 
se = translate(strel(1), [y x]);
im1 = imdilate(im1, se);


imreg = im1;
bestTransformation = [0 y x];
%%%after regiatration visualisation
if registration_plot==1
figure; visualizer(uint8(im2), im2bw(im1, graythresh(im1)), 0.4, '');title('After Registration')
s=sprintf('print -depsc %s/reg%d_im_no%d,print -djpeg %s/reg%d_im_no%d;',foldernameout,iprob,step,foldernameout,iprob,step); eval(s)

figure;imagesc(im1);colormap(gray);axis off
s=sprintf('print -depsc %s/reg%d_im1_no%d,print -djpeg %s/reg%d_im1_no%d;',foldernameout,iprob,step,foldernameout,iprob,step); eval(s)
figure;imagesc(im2);colormap(gray);axis off
s=sprintf('print -depsc %s/reg%d_im2_no%d,print -djpeg %s/reg%d_im2_no%d;',foldernameout,iprob,step,foldernameout,iprob,step); eval(s)
end


