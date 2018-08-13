function plot_matchedSpines(slice_no,slices_all,blobcentroid_all,revised_centroids,appeared_spines,lost_spines,gain_spines,foldername)

%plotting new labels on images
% figure,subplot(2,3,1),imagesc(slices_all(1).h),colormap(gray),title(num2str(1)),hold on; axis off
figure,imagesc(slices_all(1).h),colormap(gray),title(num2str(1)),hold on; title('First time point detected spines'),axis off
%We dont change labels at first time point
for k = 1:size(blobcentroid_all(1).h,1)
    text(blobcentroid_all(1).h(k,1),blobcentroid_all(1).h(k,2),num2str(k),'FontSize',8,'FontWeight','Bold','Color','y');
end
s=sprintf('print -depsc %s/first_spines_%d,print -djpeg %s/first_spines_%d;',foldername,1,foldername,1); eval(s);

ww = 10;
w = round(ww/2);

for tp = 2:slice_no %for the rest of the images through whole time sequences
% subplot(2,3,tp),imagesc(slices_all(tp).h),colormap(gray),title(num2str(tp)),hold on;axis off

im = slices_all(tp).h;
figure,imagesc(slices_all(tp).h),colormap(gray),title(['Alignment for slice number' num2str(tp)']),
hold on;axis off
for cs = 1:size(revised_centroids(tp).h,1)
    x_i = revised_centroids(tp).h(cs, 2);
    y_i = revised_centroids(tp).h(cs, 1);
    if x_i-w>=1 && y_i-w>=1 && x_i+w<=size(im,1)-1 && y_i+w<=size(im,2)-1
        if ismember(cs,appeared_spines(tp).h)
            text(revised_centroids(tp).h(cs,1),revised_centroids(tp).h(cs,2),num2str(cs),'FontSize',8, 'FontWeight','Bold','Color','c');
        else if ismember(cs,lost_spines(tp).h)
                text(revised_centroids(tp).h(cs,1),revised_centroids(tp).h(cs,2),num2str(cs),'FontSize',8, 'FontWeight','Bold','Color','r');
            else if ismember(cs,gain_spines(tp).h)
                    text(revised_centroids(tp).h(cs,1),revised_centroids(tp).h(cs,2),num2str(cs),'FontSize',8, 'FontWeight','Bold','Color','g');
                else
                    text(revised_centroids(tp).h(cs,1),revised_centroids(tp).h(cs,2),num2str(cs),'FontSize',8, 'FontWeight','Bold','Color','y');
                end
            end
        end
    end
end

s=sprintf('print -depsc %s/aligned_spines_%d,print -djpeg %s/aligned_spines_%d;',foldername,tp,foldername,tp); eval(s);
end

% Discard detected spines near the margin of the image
% ww = 40;
% w = round(ww/2);
% C1 = current_centroid;
% C2 = new_centroid;
% k = 1;
% for i = 1:size(C1,1)
%     x_i = C1(i, 2);
%     y_i = C1(i, 1);
%     if x_i-w>=1 && y_i-w>=1 && x_i+w<=size(im1,1)-1 && y_i+w<=size(im1,2)-1
%       C(k,:) = C1(i,:);
%       k = k+1;
%     end
% end
% current_centroid = C;
% clear C;
% k = 1;
% for i = 1:size(C2,1)
%     x_i = C2(i, 2);
%     y_i = C2(i, 1);
%     if x_i-w>=1 && y_i-w>=1 && x_i+w<=size(im2,1)-1 && y_i+w<=size(im2,2)-1
%       C(k,:) = C2(i,:);
%       k = k+1;
%     end
% end
% new_centroid = C;
% clear C;

