function [imTempWaterShed,imTemp,Img0,ph0,ph,waterShedLevels,r,c,c1,c2]=watershed_cv_grad_based1(imTempWaterShed,imTemp,I, filteredBinaryImage,filterBinary1,filterIndices, imMedian,numberOfComponents,maxNumberOfPixels,minNumberOfPixels,L,iprob)
%temp=zeros(size(L));

 waterShedSegmentation_all=zeros(size(I));
 waterShedSeg=zeros(size(I));
for i = 1:numberOfComponents % i=0 is the watershed, ignore
%nnz - number of nonzero elements
    if((nnz(L == i) < maxNumberOfPixels) && (nnz(L == i) > minNumberOfPixels))
            counter = i;
            
            %imTemp1 = zeros(size(imMedian));
            imTempWaterShed1 = zeros(size(imMedian));
            waterShedSegmentation = (L == i);
            imTempWaterShed1(find(waterShedSegmentation == 1)) = 255;
            tempIndices = find(imTempWaterShed1 == 255);
            waterShedLevels(i).h = imTempWaterShed1; 
            waterShedSeg=waterShedSeg+waterShedSegmentation;
            [r,c,v] = find(imTempWaterShed1);
            if(~isempty(intersect(filterIndices, tempIndices)))
                imTempWaterShed(find(waterShedSegmentation == 1)) = 255;
                waterShedSegmentation_all=waterShedSegmentation_all+waterShedSegmentation;
            end
    end
end
% beta=1.0e-6;h1=1;h2=1;maxit=100;epsilon=1;dt=0.001;alpha=0.01;
% mu=500;lambda1=1;lambda2=lambda1;Hind=2;
beta=1.0e-6;h1=1;h2=1;maxit=100;epsilon=1;dt=0.01;alpha=0.01;
mu=800;lambda1=1;lambda2=lambda1;Hind=2;
  a=waterShedSegmentation_all.* double(imMedian);
     I0= filterBinary1.* double(I); 
%        m1=mean(mean(a(find(a>0))));Img0=a;
%       Img0(find(a==0))=m1;Img0=double(Img0(:,:,1));
      Img0=double(a);
       ph0=bwdist(filteredBinaryImage); ph0=-double(ph0); %so the centers of binary image blobs are so bright
     % ph0=bwdist(waterShedSegmentation_all); ph0=-double(ph0);
   
     %[ph]=cv_grad(ph0,Img0,I0,Hind,dt,mu,alpha,lambda1,lambda2,maxit,epsilon,beta,iprob,h1,h2);

    % [ph]=cv_grad_tm(Img0,I0,ph0,beta,mu,h1,h2,maxit,epsilon,dt,lambda1,lambda2);
    
    [ph,u,c1,c2]=relax_coarsest2(Img0,ph0,beta,mu,maxit,epsilon,dt,lambda1,lambda2,2,100);
 
         imTemp(find(ph > 0)) = 255;
