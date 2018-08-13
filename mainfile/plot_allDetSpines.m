function plot_allDetSpines(slices_all,revised_centroids,foldername)

ww = 10;
w = round(ww/2);
time_point = size(revised_centroids,2);
%this is for 8th dataset which was manually segmented and labeled by Ozgur
% actual_spines = [1,2,3,5,7,8,9,10,11,12,13,16,17,18,20]; 
im = slices_all(time_point).h;
figure,imagesc(im),colormap(gray),axis off,hold on; title('All detected spines'),
for k = 1:size(revised_centroids(time_point).h,1)
    x_i = revised_centroids(time_point).h(k, 2);
    y_i = revised_centroids(time_point).h(k, 1);
    if x_i-w>=1 && y_i-w>=1 && x_i+w<=size(im,1)-1 && y_i+w<=size(im,2)-1
    %the codes in the comment are for actual spines labeled by Ozgur in 8th
    %dataset
%       if ismember(k,actual_spines)
%         text(revised_centroids(time_point).h(k,1),revised_centroids(time_point).h(k,2),num2str(k),'FontSize',12,'FontWeight','Bold','Color','g');
%       else
        text(revised_centroids(time_point).h(k,1),revised_centroids(time_point).h(k,2),num2str(k),'FontSize',8,'FontWeight','Bold','Color','y');
%       end
    end
end
s=sprintf('print -depsc %s/allDetSpines,print -djpeg %s/allDetSpines;',foldername,foldername); eval(s);
