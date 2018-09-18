clc
clear all;
close all;
clear all;
clc;

%parentdir = 'C:\desctop_NEW\work+new\Segmentation+new\Filters+SIFT+Det_Track_Seg_Bike\1data1\Ali_last';

 Im=[];
for i=1:32,     
   if i<10;
       imname = sprintf('ZSeries-neuron2time24-138_Cycle00001_CurrentSettings_Ch1_00000%d.tif',i) ;
   else
       imname = sprintf('ZSeries-neuron2time24-138_Cycle00001_CurrentSettings_Ch1_0000%d.tif',i) ;
   end 
Icube=imread(imname);
Im(:,:,i)=Icube;
%load ( sprintf('raw_%d.mat',i));

end
MIP3 = (max(Im,[],3));  % In z-direction
      MIP3 = uint8((double(MIP3) ./ 4096) .* 255); 

    imwrite(MIP3, sprintf('im_%d.tif',136));
