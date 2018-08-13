function []=Dynamic_Seg_Func(history_all,slices_all,revised_centroids,change_history,foldername)

close all
global Missing_spine_prob_Thre rad_rect
 %Missing_spine_prob_Thre=0.1;

s = sprintf('addpath .\\%s',foldername);eval(s);
tp = size(revised_centroids,2);
%this is for 5th dataset which was manually segmented and labeled by an expert
% actual_spines = [3,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,29]; 
ns = size(revised_centroids(tp).h,1);
mt_all = [];
[trackPer,lost_tp,firstAppear,allAppear] = Track_Accur(history_all,tp,ns); %how many times each spines are tracked among given number (tp) of time series(trackPer) and for how many time points a spine is lost (lost_tp).
trackPerBefore = trackPer;
%plotting tracking accuracy for each spine on a single image
% figure,imagesc(slices_all(tp).h),colormap(gray),hold on, axis off;
% for ss = 1:ns
%     format bank;
%     text(revised_centroids(tp).h(ss,1),revised_centroids(tp).h(ss,2),num2str(trackPerBefore(ss)),'FontSize',8, 'Color','r');
% end
% s=sprintf('print -depsc %s/trackAccuracy,print -djpeg %s/trackAccuracy;',foldername,foldername); eval(s);
% figure,imagesc(slices_all(tp).h),colormap(gray),hold on, axis off;
% for ss = 1:ns
%     format bank;
%     if trackPerBefore(ss)>=Missing_spine_prob_Thre
%        text(revised_centroids(tp).h(ss,1),revised_centroids(tp).h(ss,2),num2str(trackPerBefore(ss)),'FontSize',8, 'Color','r');
%     end
% end
% s=sprintf('print -depsc %s/trackAccuracy_2,print -djpeg %s/trackAccuracy_2;',foldername,foldername); eval(s);

%spine no vs time plot which shows spine apperance path.
x = 1:ns;
% figure
% for tt = 1:ns
%     plot(allAppear(tt,:),x(tt),'*','Color',rand(1,3),'LineWidth',3);
%     axis([0.05 tp+1 0 ns+1]); %axis([xmin xmax ymin ymax])
%     hold on;
% end
% grid on;
% set(gca,'yTick',0:1:(ns+1))
% xlabel('Time points')
% ylabel('Label of spine')
% s=sprintf('print -depsc %s/spinePath,print -djpeg %s/spinePath;',foldername,foldername); eval(s);

%if tracking accuracy is less than 1 and more than 0.2, we will go back to time points
%where the spines whose accuracy is in given range are not found. Then we
%will relax some parameters in detection&segmentation to make program to
%detect that spine.
missed_spines = [];
for kk = 1:size(trackPer,1)
    if (trackPer(kk)<1) && (trackPer(kk)>= Missing_spine_prob_Thre)
    missed_spines = [missed_spines,kk]; %which one we missed to detect
    end
end


for mm = 1:size(missed_spines,2)

    missed_tp = lost_tp(missed_spines(mm)).h; %at which time points the spine is missed
    all_found_tps(mm).h = setdiff(1:tp,missed_tp); %time points excluding the ones where the spine was not detected.
    dummy2 = zeros(size(missed_tp));
    for ll = 1:size(missed_tp,2)
        if ~ismember(missed_tp(ll)+1,missed_tp) && (missed_tp(ll)+1)<=tp
            dummy2(ll) = missed_tp(ll)+1;
        else if ~ismember(missed_tp(ll)-1,missed_tp)&& missed_tp(ll)>1 && (missed_tp(ll)-1)<=tp
                 dummy2(ll) = missed_tp(ll)-1;
            else if  ~ismember(missed_tp(ll)+2,missed_tp) && (missed_tp(ll)+2)<=tp
                     dummy2(ll) = missed_tp(ll)+2;
                else if ~ismember(missed_tp(ll)-2,missed_tp)&& missed_tp(ll)>2 && (missed_tp(ll)-2)<=tp
                         dummy2(ll) = missed_tp(ll)-2;
                    end
                end
            end
            
        end
         if dummy2(ll) == 0
             dummy2(ll) = all_found_tps(mm).h(1);
         end
    end

 found_tps(mm).h = dummy2;
 missed_tps(mm).h = missed_tp;

end

%After we identified which spines and at which time points are
%mis-detected, we will start dynamic segmentation algorithm here.
for nn =1:size(missed_spines,2)
    t_improve_score = 0; %tracking improve score to check whether any improvement in tracking acccuracy for missed spines
    for tt =  1:size(missed_tps(nn).h,2)
        mt = missed_tps(nn).h(tt); %time point where there is missed detection for nnth missed spine.
       ft = found_tps(nn).h(tt); %time point where we could detect that spine.
        Tfound = sprintf('T%d.mat',ft); %to get information from time points where that spine was found
        load(Tfound)
        im1 = Img;
       %%ph0 = ph; %getting contour information from time point where we found spine - as an initial point for time point where we couldnt find that spine
        if ft == 1
            cutcenter = revised_centroids(ft+1).h(missed_spines(nn),:);
            centroids = revised_centroids(ft+1).h;
        else 
            cutcenter = revised_centroids(ft).h(missed_spines(nn),:);
            centroids = revised_centroids(ft).h;
        end
    
        [ph_n,phbw,ph1]=DrawNewROI_circle(N,N,centroids, rad);%ph_n is the last level set given centroids
        ph0=ph1; 
        cutcenter = round(cutcenter);
        %%in case the spine center is close to the image border
        if size(Img,1) - cutcenter(1) < rad_rect 
             rad_rect1 = size(Img,1) - cutcenter(1);
        else if cutcenter(1) < rad_rect
                rad_rect1 = cutcenter(1);
            else
                 rad_rect1 = rad_rect;
            end
        end
         if size(Img,2) - cutcenter(2) < rad_rect 
             rad_rect2 = size(Img,2) - cutcenter(2);
         else if cutcenter(2) < rad_rect
                rad_rect2 = cutcenter(2);
             else
                 rad_rect2 = rad_rect;
             end
         end      

        rect = [cutcenter(1)-rad_rect1,cutcenter(2)-rad_rect2,rad_rect1*2,rad_rect2*2];
      %New we constract a small rectangle 2 by 2 to check if we do have a
      %spine in the 20 by 20 ROI
     % keyboard
        ph0_ROI = imcrop(ph0,rect); %selecting ROI
       % figure,imagesc(ph0_ROI),colormap(gray), axis off

        Tmissed = sprintf('T%d.mat',mt); %to work on time point where we couldnt find that spine
        load(Tmissed)
        phBefore = ph;
    
        backgroundInt = c2;
        beta=1.0e-6;h1=1;h2=1;maxit=50;epsilon=1;dt=0.001;alpha=0.01;
        mu=100^2;lambda1=1;lambda2=lambda1;Hind=2;
        maxNumberOfPixels = 6000;
       
        if mt == 1
            if size(revised_centroids(mt+1).h,1)<missed_spines(nn)
                check = firstAppear(missed_spines(nn)).h;
                sum1 = 0;
                for step = mt+1:check
                    sum1 = sum1+ change_history(step).h;
                end
                cutcenter = revised_centroids(check).h(missed_spines(nn),:)-sum1;
                cutcenterOld =cutcenter;
            else
                cutcenter = revised_centroids(mt+1).h(missed_spines(nn),:)-change_history(mt+1).h;
                cutcenterOld = cutcenter; 
            end
        else
            if size(revised_centroids(mt).h,1)<missed_spines(nn)
                check = firstAppear(missed_spines(nn)).h;
                sum1 = 0;
                for step = mt+1:check
                    sum1 = sum1+ change_history(step).h;
                end
                cutcenter = revised_centroids(check).h(missed_spines(nn),:)-sum1;
            else
                cutcenter = revised_centroids(mt).h(missed_spines(nn),:);
            end
            cutcenterOld = cutcenter;
        end
        cutcenter = round(cutcenter);
        %%in case the spine center is close to the image border
        if size(Img,1) - cutcenter(1) < rad_rect 
             rad_rect1 = size(Img,1) - cutcenter(1);
        end
        if cutcenter(1) < rad_rect
            rad_rect1 = cutcenter(1);
        end
         if size(Img,2) - cutcenter(2) < rad_rect 
             rad_rect2 = size(Img,2) - cutcenter(2);
        end
        if cutcenter(2) < rad_rect
            rad_rect2 = cutcenter(2);
        end
        %
        rect2 = [cutcenter(1)-rad_rect1,cutcenter(2)-rad_rect2,rad_rect1*2,rad_rect2*2];
    
        [Img0,L,waterShedLevels]=watershed_seg(Img,maxNumberOfPixels);
       
%         figure,imagesc(Img0),colormap(gray), axis off,hold on,plot(cutcenter(1),cutcenter(2),'rx');
%         s=sprintf('print -depsc %s/watershed_%d_%d,print -djpeg %s/watershed_%d_%d;',foldername,missed_spines(nn),mt,foldername,missed_spines(nn),mt); eval(s);
        if cutcenter(1) > size(Img0,1) || cutcenter(2) > size(Img0,2) || cutcenter(1)<=0 || cutcenter(2)<=0%which means the coordinate is out of the image borders.In this case, we skip this time point and go to next missed time point
            continue; 
        else
             
            component_no = L(cutcenter(2),cutcenter(1));
           
            if component_no == 0
               % keyboard
                continue;
            end
          
            waterShedRegion = waterShedLevels(component_no).h;
            if isempty(waterShedRegion)
               % [ph_cp,ph1_cp]=DrawNew_component_circle(n,m,cutcenter, 10);
                maxNumberOfPixels = 80000;%%%more enhancmant to get a new watershed region 
                [Img0,L,waterShedLevels]=watershed_seg(Img,maxNumberOfPixels);
                component_no = L(cutcenter(2),cutcenter(1));
                 if component_no == 0
                continue;
                end
                waterShedRegion = waterShedLevels(component_no).h;
                %%%%%show the watershade image and the missing spine centroid
%                 figure,imagesc(Img0),colormap(gray), axis off,hold on,plot(cutcenter(1),cutcenter(2),'rx');
%                 s=sprintf('print -depsc %s/watershed_W_%d_%d,print -djpeg %s/watershed_W-%d_%d;',foldername,missed_spines(nn),mt,foldername,missed_spines(nn),mt); eval(s);
            %%%new Ladi: Trying to enlarge the watershed region by threshold changes
               if isempty(waterShedRegion)
                maxNumberOfPixels = 80000;%%%more enhancmant to get a new watershed region 
                [Img0,L,waterShedLevels]=watershed_seg1(Img,maxNumberOfPixels);
                component_no = L(cutcenter(2),cutcenter(1));
                if component_no == 0
                continue;
                end
                waterShedRegion = waterShedLevels(component_no).h;
                %%%%show watershad and segmented spine
%                 figure,imagesc(Img0),colormap(gray), axis off,hold on,plot(cutcenter(1),cutcenter(2),'rx');
%                s=sprintf('print -depsc %s/watershed_%d_%d,print -djpeg %s/watershed_%d_%d;',foldername,missed_spines(nn),mt,foldername,missed_spines(nn),mt); eval(s);
                end
            end
            if isempty(waterShedRegion)
               % keyboard
                continue;
            end
        end
       
        Img0 = waterShedRegion; 
        Img0_ROI = imcrop(Img0,rect2);
        % % %         %%%%%old Bike%%%%%%%%%%%%%%%%
% %          Img0_ROI((Img0_ROI == 0)) = backgroundInt;
% %         roiImg = double(imcrop(Img,rect2));
% %         a_ROI=Img0_ROI.* double(roiImg);
% %         Img0_ROI=double(a_ROI);

%         % New%%%%%%%%%%%%%
         ph_crop=imcrop(ph,rect2);
         [m1,n1]=size(ph_crop);ph_crop_pozitive=zeros(m1,n1);
           for ph_in1=1:m1
            for ph_in2=1:n1
                if ph_crop(ph_in1,ph_in2)<=0;
                    ph_crop_pozitive(ph_in1,ph_in2)=1;
                end
            end
           end

        Img0_ROI((Img0_ROI == 0)) = backgroundInt;
        Img0_ROI((ph_crop_pozitive == 0)) = backgroundInt;
        roiImg = double(imcrop(Img,rect2));
        a_ROI=Img0_ROI.* double(roiImg);
        Img0_ROI=double(a_ROI);
  
        test_mean=sum(sum(roiImg)); 
        ts= test_mean/((rad_rect+1)^2);
        if ts<(c1+c2)/4;
         continue;
        end

        %check size of ph0_ROI and Img0_ROI cause they might be different
        %due to their center locations
        if size(Img0_ROI,1) ~= size(ph0_ROI,1) || size(Img0_ROI,2) ~= size(ph0_ROI,2)
            [n1,m1] = size(Img0_ROI);
            [n2,m2] = size(ph0_ROI);
            if n1 < n2
                rev_n = n1;
                ph0_ROI = imresize(ph0_ROI,[rev_n,m2]);
            elseif n2 <n1
                rev_n = n2;
                Img0_ROI = imresize(Img0_ROI,[rev_n,m1]);
            end
            if m1<m2
                rev_m = m1;
                ph0_ROI = imresize(ph0_ROI,[size(ph0_ROI,1),rev_m]);
            elseif m2<m1
                rev_m = m2;
                Img0_ROI = imresize(Img0_ROI,[size(Img0_ROI,1),rev_m]);
            end
        end
   %dbstop if error
   % keyboard
 [m1,n1]=size(Img0_ROI);
 if m1>3 && n1>3
    % if  mt==19 && missed_spines(nn)==9;keyboard;end%%
        [ph_new,u_new]=relax_coarsest2(Img0_ROI,ph0_ROI,beta,mu,maxit,epsilon,dt,lambda1,lambda2,1,1000);
 end
        ch = ismember(find(ph_new>0),find(Img0_ROI == 0));
        [ind,var]= find(ch);
        formatSpec = 'The %dth spine which missed at time point %d';
        str = sprintf(formatSpec,missed_spines(nn),mt);
        %%%%spine missing in a roi
        %figure;imagesc(roiImg);axis off;hold on;contour(ph_new,[0,0],'yo');colormap(gray);title(str);
        empMtx = zeros(size(Img));
        empMtx(:,:) = -1;
        if nn==2
     % keyboard
        end
        if cutcenterOld(1) <rad_rect && cutcenterOld(2) <rad_rect
           print('both...')
        elseif cutcenterOld(1) <rad_rect%%%%needsto be fixed
    %empMtx(cutcenterOld(2)-(cutcenterOld(1)-1):cutcenterOld(2)+(cutcenterOld(1)-1),1:cutcenterOld(1)+(cutcenterOld(1)-1)) = ph_new(1:(2*cutcenterOld(1)-1),1:(2*cutcenterOld(1)-1));
        elseif cutcenterOld(2) <rad_rect    
% % %     dbstop if error
% % %      diff_rad= rad_rect-cutcenterOld(2); 
% % %      empMtx(1:(2*(cutcenterOld(2)-2*diff_rad)-1),cutcenterOld(1)-(cutcenterOld(2)-1):cutcenterOld(1)+(cutcenterOld(2)+1)) = ph_new(1:(2*(cutcenterOld(2)-2*diff_rad)-1),1:(2*cutcenterOld(2)+2*diff_rad-1));
   % empMtx(1:cutcenterOld(2)+(cutcenterOld(2)-1),cutcenterOld(1)-(cutcenterOld(2)-1):cutcenterOld(1)+(cutcenterOld(2)-1)) = ph_new(1:(2*cutcenterOld(2)-1),1:(2*cutcenterOld(2)-1));
           
        else
            empMtx(cutcenterOld(2)-(size(ph_new,1)-((size(ph_new,1)+1)/2)):cutcenterOld(2)+(size(ph_new,1)-((size(ph_new,1)+1)/2)),cutcenterOld(1)-(size(ph_new,2)-((size(ph_new,2)+1)/2)):cutcenterOld(1)+(size(ph_new,2)-((size(ph_new,2)+1)/2))) = ph_new;
        end
            dbstop if error
        reseg_spine(nn).time(mt).h = empMtx;
        formatSpec = 'Spine(s) missed at time point %d';
       
        str = sprintf(formatSpec,mt);
        if nn == 1
            figure,imagesc(Img),axis off,hold on,contour(phBefore,[0,0],'yo'),colormap(gray);title(str);
            hold on,
            contour(empMtx,[0,0],'ro');
            s=sprintf('print -depsc %s/redetectWhole_%d,print -djpeg %s/redetectWhole_%d;',foldername,mt,foldername,mt); eval(s);
        else
            if  mt_all ~= 0 
                if ismember(mt,mt_all(1,:))
                    id = find(mt_all(1,:) == mt);
                    figure,imagesc(Img),axis off,hold on,contour(phBefore,[0,0],'yo'),colormap(gray);
                    contour(empMtx,[0,0],'ro'),hold on,
                    for ii = 1:size(id,2)
                        contour(reseg_spine(mt_all(2,id(ii))).time(mt).h,[0,0],'ro');
                        hold on
                    end
                    colormap(gray); title(str);
                    s=sprintf('print -depsc %s/redetectWhole_%d,print -djpeg %s/redetectWhole_%d;',foldername,mt,foldername,mt); eval(s);
                else
                    figure,imagesc(Img),axis off,hold on,contour(phBefore,[0,0],'yo'),colormap(gray);title(str);
                    hold on,
                    contour(empMtx,[0,0],'ro');
                    s=sprintf('print -depsc %s/redetectWhole_%d,print -djpeg %s/redetectWhole_%d;',foldername,mt,foldername,mt); eval(s);
            
                end
            else
                figure,imagesc(Img),axis off,hold on,contour(phBefore,[0,0],'yo'),colormap(gray);title(str);
                hold on,
                contour(empMtx,[0,0],'ro');
                s=sprintf('print -depsc %s/redetectWhole_%d,print -djpeg %s/redetectWhole_%d;',foldername,mt,foldername,mt); eval(s);
            end
        end

         mt_all = [mt_all,vertcat(mt,nn)];
         if size(ind,1) ~= size(find(ph_new>0),1) %which means we could segment the spine now.
            t_improve_score = t_improve_score+1;
        elseif size(find(empMtx>0),1) == size(find(ph_new>0),1) %which means we could segment the spine now.
            t_improve_score = t_improve_score+1;
        else
            %we need shape prior information to improve segmentation
        end
    end
    trackPer(missed_spines(nn)) = trackPer(missed_spines(nn))+t_improve_score/tp;
    close all
end

% figure;imagesc(slices_all(tp).h);colormap(gray);hold on; title('Tracking accuracy after re-segmentation');axis off;
% for ss = 1:ns
%     format bank;
%     text(revised_centroids(tp).h(ss,1),revised_centroids(tp).h(ss,2),num2str(trackPer(ss)),'FontSize',8, 'Color','r');
% end
% s=sprintf('print -depsc %s/trackAccuracyImproved,print -djpeg %s/trackAccuracyImproved;',foldername,foldername); eval(s);

figure;imagesc(slices_all(tp).h);colormap(gray);hold on; title('Tracking accuracy after re-segmentation');axis off;
for ss = 1:ns
    format bank;
    if trackPer(ss)>=Missing_spine_prob_Thre
       text(revised_centroids(tp).h(ss,1),revised_centroids(tp).h(ss,2),num2str(trackPer(ss)),'FontSize',8, 'Color','r');
    end
end
s=sprintf('print -depsc %s/trackAccuracyImproved_2,print -djpeg %s/trackAccuracyImproved_2;',foldername,foldername); eval(s);

