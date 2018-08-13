%% Detection and segmentation ended. In below, tracking and association start.
function [history_all,slices_all,revised_centroids,change_history] = Tracking_Func_Features(foldername,iprob)

s = sprintf('addpath .\\%s',foldername);eval(s);
sliceno = iprob;
for ka=1:sliceno
T = sprintf('T%d.mat',ka); 

load(T)
centroids = MyCentroids(numberOfBlobs,blobMeasurements);
blobcentroid_all(ka).h= centroids; %saves centers' coordinates of current image
blobfeatures_all(ka).h = features;
slices_all(ka).h = Img;
descriptors(ka).h = desc;
locations(ka).h = locs;
groups(ka).h = group;
end

current_state=blobfeatures_all(1).h(:,1:5);
new_state=blobfeatures_all(2).h(:,1:5);
history_all = [];

for step = 2:(sliceno)
    
im1 = slices_all(step-1).h;
im2 = slices_all(step).h;
desc1 = descriptors(step-1).h;
locs1 = locations(step-1).h;
desc2 = descriptors(step).h;
locs2 = locations(step).h;
[change,im_reg] = landmarkBasedRegistration(im1,im2,desc1,locs1,desc2,locs2,foldername,step,iprob);
x_shift = change(:,size(change,2));
y_shift = change(:,(size(change,2)-1));
change_history(step).h = [x_shift y_shift];

current_state_new = current_state;
for i = 1:size(current_state,1)
current_state_new(i,1:2) = current_state(i,1:2)+[x_shift y_shift];
end

current_state = current_state_new;
current_state_new = [];

[m1,n1]=size(current_state);
[m2,n2]=size(new_state);
dist_history = zeros(m2,m1);
dist_history_10 = zeros(m2,m1);
new_spine_id = [];
lost_spine_id = [];

all_points_vect = zeros(m1,7);
%%distance measurement  
for ia=1:m1 %t time frame
    r0=current_state(ia,:); %takes iath spine's center coordinates of first image
for ib=1:m2 %t+1 time frame
      r1=new_state(ib,:); %takes ibth spine's center coordinates of second image
      dist=sqrt((r1(:,1)-r0(:,1))^2+(r1(:,2)-r0(:,2))^2+(r1(:,4)-r0(:,4))^2+(r1(:,5)-r0(:,5))^2); %finds euclidean distance between spine centers of new image and previous image
      dist_vect(1,ia)=dist; %stores distances of each spines
      dist_history(ib,ia) = dist;
      if dist<20
          %distance threshold
        dist_prev = dist_history_10(:,ia);
        if ib >1
            if find(dist_prev == 1) ~= 0 %which means there was another time point which satisfied threshold before
                prev_ones = find(dist_prev == 1);
                if dist_history(find(dist_prev == 1),ia) < dist %if previous point has smaller distance which satisfies the threshold
                    dist_vect(2,ia)=0;  %it means - this spine is common in each image
                else
                    dist_vect(2,ia)=1; % it means - this spine is new or missing
                    all_points_vect(ia,1:7)=[r0(:,1),r0(:,2),r1(:,1),r1(:,2),ia,ib,dist]; %it creates a vector storing spine's coordinates/common spine numbers/ distances between common spines
                end
            else
                dist_vect(2,ia)=1;  %it means - this spine is common in each image
                all_points_vect(ia,1:7)=[r0(:,1),r0(:,2),r1(:,1),r1(:,2),ia,ib,dist]; %it creates a vector storing spine's coordinates/common spine numbers/ distances between common spines
            end
        else
            dist_vect(2,ia)=1;  %it means - this spine is common in each image
            all_points_vect(ia,1:7)=[r0(:,1),r0(:,2),r1(:,1),r1(:,2),ia,ib,dist]; %it creates a vector storing spine's coordinates/common spine numbers/ distances between common spines
        end
      else
        dist_vect(2,ia)=0; % it means - this spine is new or missing
      end 
      dist_history_10(ib,ia) = dist_vect(2,ia);          
end
  
end

% to erase duplication due to having more than one points which fits to
% distance threhold (we choose the one which has smaller distance than the
% other)
for ic = 1:size(all_points_vect,1)
    for id = 1: size(all_points_vect,1)
        if (ic ~= id)
           if (all_points_vect(ic,6) == all_points_vect(id,6))
             if (all_points_vect(ic,7) > all_points_vect(id,7))
                 all_points_vect(ic,:) = zeros(1,7);
                 dist_history_10(:,ic) = zeros(size(dist_history_10,1),1);
             else
                 all_points_vect(id,:) = zeros(1,7);
                 dist_history_10(:,id) = zeros(size(dist_history_10,1),1);
             end
           end
        end
    end
end

for id1 = 1:size(dist_history_10,1)
    if(dist_history_10(id1,:) == zeros(1,size(dist_history_10,2)))
         new_spine_id = [new_spine_id,id1];
    end
end

for id2 = 1:size(dist_history_10,2)
    if(dist_history_10(:,id2) == zeros(size(dist_history_10,1),1))
         lost_spine_id = [lost_spine_id,id2];
    end
end

%Saving outputs
dist_for_all_images(step-1).h=dist_history;
all_points_tracked(step-1).h=all_points_vect;
new_spines(step).h = new_spine_id;

%%We are done with tracking and matching spines with respect to their
%%centroids location. Next step is to associate matched spines.

%initialization of parameters
id_curr = []; %which shows if detected spines in t frame corresponds to the one in t+1 frame
id_new = []; %which shows if detected spines in t+1 frame corresponds to the one in t frame
changed_labels = []; %this will be used to revise ordering of centroid coordinates for next time point
lost_spine_id2 = [];
appearspinesid = [];
gain_spine_id = [];
gain_spine_label = [];
myvect = all_points_vect(:,5:6);
myvect2 = zeros(max(max(myvect)),2);

for step2 = 1:size(myvect,1)
%steady and lost spines matching and relabelling
    if myvect(step2,1) == myvect(step2,2) && (myvect(step2,2) ~= 0)
       id_curr = [id_curr;1];
       id_new = [id_new;1];
       myvect2(step2,1) = myvect(step2,2);
       myvect2(step2,2) = myvect(step2,1);
       changed_labels = [changed_labels;horzcat(myvect2(step2,1), myvect2(step2,2))];
    else if  myvect(step2,2) == 0 %there is a lost spine or reappeared spine at t+1 frame 
             if step >2
                 check_prev_lost = lost_spines(step-1).h;
                 if  ismember(step2,lost_spine_id) %which means this one is lost spine
                     id_curr = [id_curr;1];
                     id_new = [id_new;-1];
                     myvect2(step2,1) = 0; 
                     myvect2(step2,2) = step2;    
                 end
             else %which means this one is lost spine
                 id_curr = [id_curr;1];
                 id_new = [id_new;-1];
                 myvect2(step2,1) = 0; 
                 myvect2(step2,2) = step2; 
                 lost_spine_id2 = [lost_spine_id2,step2];
             end
        else 
             id_curr = [id_curr;1];
             id_new = [id_new;1];
             myvect2(step2,1) = myvect(step2,2);
             myvect2(step2,2) = myvect(step2,1); 
             changed_labels = [changed_labels;horzcat(myvect2(step2,1), myvect2(step2,2))];
        end
    end
end
history = [id_curr,id_new]; 
%checking previous "spine appearence history list" after first step
if step >2
prev_history_check = history_all(1:size(history,1),step-1);
for hs = 1:size(prev_history_check,1)
    if prev_history_check(hs,:) == -1 && history(hs,2) == 1
        id_curr = [id_curr;-1];
        id_new = [id_new;1];
        appearspinesid = [appearspinesid;hs];
    end
end
end 

checkvec = 1:max(size(new_state,1));
nocorr = find(ismember(checkvec',myvect2(:,1)) == 0)'; %to check if there is missing number in matched spines at t+1 time point
count = 0;

if nocorr ~= 0 %which means I have new spines at t+1 frame
    num_new = size(nocorr,2); %number of new spines at t+1 frame (excluding appeared spines)
    updatehistory =zeros(size(myvect,1)+num_new,2);
    updatehistory(1:size(history,1),:) = history;
    updatehistory(size(updatehistory,1)-num_new+1:end,:) = repmat([-1,1],num_new,1);
    for p = 1:size(nocorr,2)
        id = nocorr(:,p);
        gain_spine_id = [gain_spine_id,id];
        count = count +1; 
        gain_spine_label = [gain_spine_label, size(updatehistory,1)-num_new+count];
    end
    history_nl(step).h = updatehistory;
    %revise changed_labels and myvect2
    gain_spine_id = sort(gain_spine_id,'descend');
    for gc = 1:size(gain_spine_id,2)
        changed_labels = [changed_labels;horzcat(gain_spine_id(gc),gain_spine_label(gc))];
    end
else
    history_nl(step).h = history;
end
revised_labels(step).h = changed_labels;
lost_spines(step).h = lost_spine_id;
appeared_spines(step).h = appearspinesid; 
gain_spines(step).h = gain_spine_label;

%We created a matrix called history_all which we can see at which time point
%each spine exists, disappears and reappears. 
if step == 2
    history_all = history_nl(step).h;
else
    if size(history_nl(step).h,1) > size(history_nl(step-1).h,1)
        history_all2 = zeros(size(history_nl(step).h,1),step);
        differ = size(history_nl(step).h,1) - size(history_all,1);
        history_all2(1:size(history_all,1),:) = horzcat(history_all,history_nl(step).h(1:size(history_all,1),2));
        history_all2(size(history_all2,1)-differ+1:end,:) = horzcat(repmat(-1,differ,step-1),history_nl(step).h(size(history_all2,1)-differ+1:end,2));
    else if size(history_nl(step).h,1) < size(history_nl(step-1).h,1)
            history_all2 = zeros(size(history_all,1),step);
            differ = size(history_all,1) - size(history_nl(step).h,1);
            history_all2(1:size(history_nl(step).h,1),:) = horzcat(history_all(1:size(history_nl(step).h,1),:),history_nl(step).h(:,2));
            history_all2(size(history_all,1)-differ+1:end,:) = horzcat(history_all(size(history_all,1)-differ+1:end,:),repmat(-1,differ,step-1));
        else
            history_all2 = horzcat(history_all,history_nl(step).h(:,2));
        end
        
    end
    history_all = history_all2;
end

%Defining spines' centroids for next step alignment & matching process
if step == sliceno
    current_state = blobfeatures_all(step).h(:,1:5);
else
    new_state = blobfeatures_all(step+1).h(:,1:5);
    current_state = blobfeatures_all(step).h(:,1:5);
end

%changing the locations of centroids based on relabelling results
new_cc = zeros(size(history_nl(step).h,1),5); %new current centroid
checkvec2 = 1:max(changed_labels(:,1));
nocorr2 = find(ismember(checkvec2',changed_labels(:,1)) == 0)'; 
for c =1:size(changed_labels,1)
    new_cc(changed_labels(c,2),:) = current_state(changed_labels(c,1),:);
end
if nocorr2 ~= 0 %this step is to capture new or reappeared spines
    [lia,lib] = ismember(nocorr2,myvect2(:,1));
    new_cc(myvect2(lib,2),:) = current_state(myvect2(lib,1),:);
end
if step>2
    if lost_spines(step).h ~=0
        for ad = 1:size(lost_spines(step).h,2)
            %checking where this lost spine was appeared before it got lost
            %until this time point
            count2 = 0;
            while step-count2>1
               if history_all(lost_spines(step).h(ad),step-count2-1) == 1
                 count2 = count2+1;
                 break
               else
                 count2 = count2+1;
                 continue
               end
              
            end
            if step-count2 == 1 %if that spine was appeared only at first time point
               new_cc(lost_spines(step).h(ad),:) = blobfeatures_all(step-count2).h(lost_spines(step).h(ad),1:5); 
            else
               new_cc(lost_spines(step).h(ad),:) = blobfeatures_all(step-count2).h(revised_labels(step-count2).h(find(revised_labels(step-count2).h(:,2) == lost_spines(step).h(ad)),1),1:5);
            end
            if count2 >1
                sum = 0;
                for step3 = step-count2+1:step
                    sum = sum+ change_history(step3).h;
                end
                new_cc(lost_spines(step).h(ad),:) = [new_cc(lost_spines(step).h(ad),1:2)+sum,new_cc(lost_spines(step).h(ad),3:5)];
            else if count2 == 1 
%                     && step == sliceno
                    new_cc(lost_spines(step).h(ad),:) = [new_cc(lost_spines(step).h(ad),1:2)+change_history(step).h,new_cc(lost_spines(step).h(ad),3:5)];
                end
            end
        end
    end
else
    if lost_spines(step).h ~=0
        for ad = 1:size(lost_spines(step).h,2)
            new_cc(lost_spines(step).h(ad),:) = [blobfeatures_all(step-1).h(lost_spines(step).h(ad),1:2)+change_history(step).h, blobfeatures_all(step-1).h(lost_spines(step).h(ad),3:5)];
        end
    end
end

current_state = new_cc;
revised_centroids(step).h = current_state(:,1:2);


clear all_points_vect dist_vect 

end
plot_matchedSpines(sliceno,slices_all,blobcentroid_all,revised_centroids,appeared_spines,lost_spines,gain_spines,foldername) %this function displays results on time series
plot_allDetSpines(slices_all,revised_centroids,foldername);  %this function displays all detected spines on one final image



