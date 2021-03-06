clear all; 
close all; 
clc;     
  
% addpath   ..\1data1/2014-01-31-s1 ;  dataset=1; %slice_no=29
addpath  ..\1data1/9sept2011;       dataset=2;%PixelType=3;slice_no=53
%addpath  ..\1data1/13April2012;     dataset=3;%PixelType=3;slice_no=41
%addpath  ..\1data1/Data4;          dataset=4;%slice_no=52
%addpath  ..\1data1/Data5;           dataset=5;%slice_no=10
%addpath  ..\1data1/26Feb2011;       dataset=6;%PixelType=3;,slice_no=41
%addpath  ..\1data1/29Jan2011;       dataset=7;%slice_no=53
%addpath  ..\1data1/11June2011;      dataset=8;%slice_no=56
%addpath  ..\New-data-set-8TrackingData4Times/MIP_FOR_ALL;      dataset=9;%slice_no=4 wich corespond to time 1 10 20 and 30
% addpath  ..\1data1/Jochen_data;           dataset=10;%slice_no=8
%addpath  ..\1data1/Raw;   dataset=11;%slice_no=24
%addpath  ..\New-data-set-9TrackingDataManyTimes/MIP_FOR_ALL;      dataset=12;%slice_no=4 wich corespond to time 1 10 20 and 30
%addpath ..\1data1/14April2012;     dataset=13;%slice_no=41
%addpath ..\1data1/Ali_new;     dataset=14;%PixelType=1;slice_no=6
% addpath ..\1data1/June302012;     dataset=15; 
%addpath ..\1data1/Jul162011;     dataset=16; 
%addpath ..\1data1/Jul162011_each3times;     dataset=17; 
%addpath ..\1data1/June302012_eah3times;     dataset=18; 
%addpath ..\1data1/June29th2012_each3times;     dataset=19; 
%addpath ..\1data1/Oct8th2011_each3times;     dataset=20; 
  
%%%% Important Note for maxNumberOfPixels = 6000 or 3000
%%%%for 9sept2011;   medfilt_size=9;Dist_Thresh=30; rad=5;rad2=60;medfilt_size=9; maxNumberOfPixels = 3000; minNumberOfPixels = 20;
%%%13April2012; medfilt_size=9;Dist_Thresh=30; rad=5;rad2=60;medfilt_size=9; maxNumberOfPixels = 3000; minNumberOfPixels = 20;
%%%%%%%%%%%2014-01-31-s1;Data4;Raw;29Jan2011;New-data-set-8TrackingData4Times; we use dotEnh=1;alpha=5;Missing_spine_prob_Thre=0.2; Dist_Thresh=30;medfilt_size=9; maxNumberOfPixels = 3000; minNumberOfPixels = 20;
%%%%for 26Feb2011; medfilt_size=9;Dist_Thresh=30; rad=5;rad2=60;medfilt_size=9; maxNumberOfPixels = 6000; minNumberOfPixels = 20;
%%%%%%%%11June2011;Raw;29Jan2011;New-data-set-8TrackingData4Times;  we use dotEnh=1;alpha=5;Missing_spine_prob_Thre=0.2; Dist_Thresh=30;medfilt_size=9; maxNumberOfPixels = 6000; minNumberOfPixels = 20;


 load descAll.mat;
load labelsAll.mat;
%  load descAll_5training.mat;
%  load labelsAll_5training.mat;
%   load desc9Sept_im2_15_30_54.mat;%%%%
%   load labels9Sept_im2_15_30_54.mat;
%  load desc9Sept_im2_15_30_54_dotEnh.mat;
%   load labels9Sept_im2_15_30_54_dotEnh.mat;
 global rad rad2 dmin dmax N_m Missing_spine_prob_Thre Dist_Thresh Loc_Registration_swich_on rad_rect registration_plot  plot_matchedSp
 

 %%%%%filter multilevep parameters
 registration_plot=0; %to check the registration proces swich registration_plot to 1
 plot_matchedSp=1; %to check the tracking process swich registration_plot to 1
 
dmin=4;dmax=32;N_m=8;  %much proper for iprob 17, 26 and 27
 rad=10;rad2=40;%from iprob 12 smaller rad=25, for iprob12 epsilon=0.01  
%  dn = dateasstring(1);
det_postdet = 'If only detection is performed press 1 else if redetection will be applied press 2:';
decD = input(det_postdet); 
prompt = 'Please enter number of time points:';
slice_no = input(prompt);    
%%%%%%%%%%%%%%%%%%
Missing_spine_prob_Thre=0.2; 
Loc_Registration_swich_on=0;%if 1 we have some b-spline registration       
%creating a folder to save output plots.
dn = dateasstring(1);
prompt = 'Please write output foldername:';
TestID = input(prompt,'s');
dotEnh=1; %  If  dotEnh=0 or alpha=0 we do not include dot enhancment 
alp = 'Please enter $\alpha$ the enhancement parameter ( 3<0 \alpha <=5 recommended ):';
alpha = input(alp);
 
if isempty(TestID)
    TestID = 'Output'; %if user doesnt enter a name, default name is given.
end
foldername = sprintf('%s_%s',TestID,dn);
mkdir(foldername) 
type_px = 'Choose the Pixel Size of the taken images: \n Type= 1 ->  Pixel Size: 67.58x67.58 micrometer (16 bit)\n  Type=2 -> Pixel Size: 31.38x31.38 micrometer (8 bit) \n Type=3 -> Pixel Size: 19.30x19.30 micrometer (12 bit in 16 bit format): ';
PixelType = input(type_px);
if PixelType==3
    medfilt_size=9;
    Dist_Thresh=30;%distance threshold if dist<35 %distance threshold = 23 for first dataset. 35 is for 7th dataset
    rad=7;rad2=55;
    rad_rect = 22;%%%%%crop rectangle radius
elseif PixelType==2
    medfilt_size=7;
    Dist_Thresh=15;
    rad=5;rad2=40;
    rad_rect = 15;%%%%%rectangle radius
else PixelType==1
    medfilt_size=5;
    Dist_Thresh=10;
    rad=2;rad2=10;
    rad_rect = 10;%%%%%rectangle radius
end


%%%% Important Note for maxNumberOfPixels = 6000 or 3000
%%%%for 9sept2011;13April2012;2014-01-31-s1;Data4;Raw;29Jan2011;New-data-set-8TrackingData4Times; we use dotEnh=1;alpha=5;Missing_spine_prob_Thre=0.2; Dist_Thresh=30;medfilt_size=9; maxNumberOfPixels = 3000; minNumberOfPixels = 20;
%%%%for 26Feb2011;11June2011;Raw;29Jan2011;New-data-set-8TrackingData4Times;  we use dotEnh=1;alpha=5;Missing_spine_prob_Thre=0.2; Dist_Thresh=30;medfilt_size=9; maxNumberOfPixels = 6000; minNumberOfPixels = 20;

 maxNumberOfPixels = 3000;%  
minNumberOfPixels = 20;   
% svmStruct = svmtrain(descAll, labelsAll);
options.MaxIter = 100000;
svmStruct = svmtrain(descAll,labelsAll, 'Options', options); 
sumAreaVec=[];sumArea2Vec=[];
for iprob =1:slice_no
    
if dataset==2;imname = sprintf('image_%d.png',iprob+1) ;
elseif dataset==3; imname = sprintf('image_%d.png',iprob) ;
elseif dataset==4; if iprob<6;imname = sprintf('ZSeries-neuron1dhpg%d-0%d.png',iprob,iprob+41);
        else imname = sprintf('ZSeries-neuron1time%d-0%d.png',iprob,iprob+41);end;
elseif dataset==5; if iprob<10;imname = sprintf('d1-00%d_Cycle00001_Ch1_MIP.tif',iprob) ;else imname = sprintf('d1-0%d_Cycle00001_Ch1_MIP.tif',iprob) ;end
elseif dataset==6; imname = sprintf('im_%d.tif',iprob) ;
 elseif dataset==7; imname = sprintf('im_%d.tif',iprob) ;
 elseif dataset==8; imname = sprintf('im_%d.tif',iprob) ;
  elseif dataset==9; imname = sprintf('d1-00%d.png',iprob) ;
  elseif dataset==10; imname = sprintf('imagesJ_%d.png',iprob) ;
  elseif dataset==11; imname = sprintf('im_%d.tif',iprob) ;
elseif dataset==14;imname = sprintf('img%d.jpg',iprob-1) ;
    elseif dataset==15; imname = sprintf('im_%d.tif',iprob+107) ;
         elseif dataset==16; imname = sprintf('im_%d.tif',iprob) ;
              elseif dataset==17; imname = sprintf('im_%d.tif',iprob) ;
                  elseif dataset==18; imname = sprintf('im_%d.tif',iprob+107) ;
                  elseif dataset==19; imname = sprintf('MIP_%d.png',iprob) ;
                     elseif dataset==20; imname = sprintf('MIP_%d.png',iprob) ;
      elseif dataset==13; imname = sprintf('image_%d.png',iprob+11) ;
   elseif dataset==12;if iprob<10;imname = sprintf('d1-00%d.jpg',iprob);else imname = sprintf('d1-0%d.jpg',iprob);end
elseif dataset==1;if iprob<10;imname = sprintf('d1-00%d.jpg',iprob);else imname = sprintf('d1-0%d.jpg',iprob);end
end  

I=imread(imname);I=I(:,:,1); [n,m]=size(I);N=n;


imMedian = medfilt2(I, [medfilt_size medfilt_size]);
Img=imMedian;

if dotEnh==1
    [zdot,zline]=dot_line_2D_enh(Img,foldername,iprob);
   a=zdot+zline;a=uint8(a);Img=Img+alpha*a;%%%Data 5 unit16
   %Img=Img*255.\max(Img(:))+a*255.\max(a(:));
end
% figure; imagesc(Img);colormap(gray) ;axis off; %plots and saves original image
% s=sprintf('print -depsc %s/Original_%d_N%d,print -djpeg %s/Original_%d_N%d;',foldername,iprob,N,foldername,iprob,N); eval(s)

%finding SIFT key features and descriptors.
[image desc locs] = sift(Img); 
 group = svmclassify(svmStruct, desc);
 descr(iprob).h = desc;
 locas(iprob).h = locs;
 spinelocs(iprob).h = [locs(find(group==1),2),locs(find(group==1),1)];
 
[filteredBinaryImage,filBinImg1]=showkeys(Img, locs, group,iprob,foldername,N);

imEq = adapthisteq(Img); %Contrast-limited adaptive histogram equalization
%region based method that connects region with same intensity
imEMaxima = imextendedmax(Img,5);%computes the extended-maxima transform, which is the regional maxima of the H-maxima transform
imEMaxima = imfill(imEMaxima, 'holes');%fills holes in the binary image BW. 

imBW = im2bw(imEq, graythresh(imEq));%Convert an image to a binary image, based on threshold

% roi selection 
imBW = imdilate(imBW | imEMaxima, ones(3,3));
% figure; imshow(imBW) 

% find boundaries
imBW = imfill(imBW,'holes');

% watersheds (segmentation phase 1)
imEqComp = imcomplement(imEq);%%%Complement image.
% make background and imEMaxima ==> only maxima; modifies the intensity imEqComp using morphological reconstruction so it only has regional minima wherever imBW is nonzero.
imEqCompModified = imimposemin(imEqComp, ~imBW | imEMaxima);%%%%%Impose minima.
% apply watershed transform.
L = watershed(imEqCompModified); %
%figure,imagesc(L)
% maxNumberOfPixels = 6000;%6000 better for first dataset.8000 better for 2nd dataset,
% minNumberOfPixels = 20;  

% graph-based (segmentation phase 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberOfComponents = max(max(L));
imMedianTemp = Img;
   
imTemp = zeros(size(Img));
imTempWaterShed = zeros(size(Img));
filterIndices = find(filteredBinaryImage == 1);

[imTempWaterShed,imTemp,Img_cut,ph_cut,ph,waterShedLevels,r,c,c1,c2]=watershed_cv_grad_based1(imTempWaterShed,imTemp,I,filteredBinaryImage,filBinImg1,filterIndices, Img,numberOfComponents,maxNumberOfPixels,minNumberOfPixels,L,iprob);

%figure,imagesc(Img_cut),colormap(gray),hold on,contour(ph_cut,[0,0],'ro')
   
imTempWaterShed = im2bw(imTempWaterShed, graythresh(imTempWaterShed));
  
imTemp = im2bw(imTemp, graythresh(imTemp));
imRemoved = bwmorph(imTemp, 'remove'); 
imMedian(find(imRemoved == 1)) = 255;
figure; imagesc(imMedian);colormap(gray),hold on, contour(imRemoved,[1 1], 'yo');axis off
s=sprintf('print -depsc %s/imMedian_%d_Ny%d,print -djpeg %s/imMedian_%d_Ny%d;',foldername,iprob,N,foldername,iprob,N); eval(s)

%%%%%%%%%%%%%%%%    

labeledImage = bwlabel(imTemp, 8);  % Label each blob so we can make measurements of it
%%%%%this figure shows the new added spines
figure; imagesc(imMedian);colormap(gray), axis off
hold on;
boundaries = bwboundaries(imTemp);	
numberOfBoundaries = size(boundaries);
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'r', 'LineWidth', 2);title('Labeled Spines')
end

blobMeasurements = regionprops(labeledImage, imMedian, 'all');   
numberOfBlobs = size(blobMeasurements, 1);
fontSize = 10;	% Used to control size of "blob number" labels put atop the image.
labelShiftX = -7;	% Used to align the labels in the centers of the coins.
for k = 1 : numberOfBlobs % Loop through all blobs.
	blobCentroid = blobMeasurements(k).Centroid; % Get centroid.
    blobRect(k,:) = blobMeasurements(k).BoundingBox;
    roundness = (4*pi*blobMeasurements(k).Area)/blobMeasurements(k).Perimeter;
    features(k,:)  = double([blobCentroid,blobMeasurements(k).Area,blobMeasurements(k).Perimeter,roundness,blobMeasurements(k).Orientation, blobMeasurements(k).Eccentricity]);
	text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(k), 'FontSize', fontSize, 'FontWeight', 'Bold', 'Color','y');
end

s=sprintf('print -depsc %s/imMedianLab_%d_Ny%d,print -djpeg %s/imMedianLab_%d_Ny%d;',foldername,iprob,N,foldername,iprob,N); eval(s)
      
s  = regionprops(labeledImage, 'centroid');
centroids = cat(1, s.Centroid);
[sumArea2,sumArea]=labeledImageArea(imMedian,labeledImage,centroids,foldername,iprob,N);
sumArea2Vec=[sumArea2Vec;sumArea2];
sumAreaVec=[sumAreaVec;sumArea];
save(sprintf('%s/T%d.mat',foldername,iprob))

clear Img  t features blobRect
close all

end

%Tracking starts - distance metric matching which compares spine centroids
%  foldername ='out_20180808'; iprob=2; %iprob is the no of slices you used
%  for the experiment 
        
[history_all,slices_all,revised_centroids,change_history] = Tracking_Func_Dist(foldername,iprob);
s=sprintf('save %s/Trac.mat history_all slices_all revised_centroids change_history foldername', foldername);eval(s)      
%Tracking starts - distance metric matching which compares spine centroids,
%area,perimeter and roundness (4*pi*area/perimeter)
% [history_all,slices_all,revised_centroids,blobrect_all,change_history] = Tracking_Func_Feature(foldername,iprob);
if decD==1
    break   
else
%Dynamic segmentation starts
Dynamic_Seg_Func1(history_all,slices_all,revised_centroids,change_history,foldername);
end 
