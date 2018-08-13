function  [zDotEach,terminating_pts1,filteredImage,filBinImg1]=dot_line_2D_rgb1(imMedian,subfoldername,iprob,rad, bWantBright)

global dmin dmax N_m
%[n,m]=size(imMedian);

% Img = imresize(double(origImg),[n,n]);%
Img = double(imMedian(:,:,1));
figure;imagesc(Img);colormap(gray);%title('Smooth image');
s=sprintf('print -depsc %s/SmoothImg_%d_Nm%d,print -djpeg %s/SmoothImg_%d_Nm%d;',subfoldername,iprob,N_m,subfoldername,iprob,N_m); eval(s)
figure;imagesc(Img);colormap(gray);axis off%title('Smooth image');
s=sprintf('print -depsc %s/nSmoothImg_%d_Nm%d,print -djpeg %s/nSmoothImg_%d_Nm%d;',subfoldername,iprob,N_m,subfoldername,iprob,N_m); eval(s)

brightFlag = 1; %%% highlight bright object
if(nargin>4)
    brightFlag = bWantBright;
end

%% loop of multi-scale

r = nthroot(dmax/dmin, N_m);
sigma = dmin/4;

zdotall = zeros(size(Img, 1), size(Img, 2))-realmin;
zlineall = zeros(size(Img, 1), size(Img, 2))-realmin;

zDotEach = zeros(N_m, size(Img, 1), size(Img, 2));
zLineEach = zeros(N_m, size(Img, 1), size(Img, 2));

rgbDebug=[]; 
for i=1:N_m
    %sigma = sigma*r^(i-1);
    sigma = sigma*r;
    G = fspecial('gaussian', ceil(4*sigma), sigma);
    smoothedImg = imfilter(Img, G, 'conv');
 

    %% calculate 1st and second derivatives.
    [fx,fy]      =   gradient(smoothedImg);
    [fxx,fxy]   =   gradient(fx);
    [fyx,fyy]   =   gradient(fy);
    
    
    %% calculate K, H, F.
    % gray level;
    K = fxx(:,:,1);
    F = fxy(:,:,1);
    H = fyy(:,:,1);

    for j=2:size(Img, 3)
    %rgb
    %if(rgbFlag)
        K = K+ fxx(:,:,j);
        F = F+ fxy(:,:,j);
        H = H+ fyy(:,:,j);
    %end
    end
    
    %% calculate eigenvalues
    P = (K+H)/2;
    %Q = sqrt(K.*H-F.*F); %element-wise multiplication.

    Q = sqrt(P.^2 - K.*H + F.*F);
    
    lamda1 = P + Q;
    lamda2 = P - Q;

    rgbDebug = [rgbDebug; max(max(K)) max(max(H)) max(max(F)) max(max(P)) max(max(Q)) max(max(lamda1)) max(max(lamda2))];
    
    clear fx fy fxx fyy fxy fyx K F H P Q
    
    [row, col] = find(abs(lamda1)<abs(lamda2));
    
    if(size(row, 1)>0)
        tmpLamda = lamda1;
        index = (col-1) * size(Img, 1) + row; % see matrix indexing
        lamda1(index) = lamda2(index);
        lamda2(index) = tmpLamda(index);        
    end

    %% compute dot
    zdot = zeros(size(Img, 1), size(Img, 2));
    
    
    if(brightFlag == 1)
        [r1, c1] = find(lamda1<0);
        [r2, c2] = find(lamda2<0);  
    else
        [r1, c1] = find(lamda1>0);
        [r2, c2] = find(lamda2>0);  
    end
    
    interRC = intersect([r1 c1], [r2 c2], 'rows');
    clear r1 c1 r2 c2
    index = (interRC(:,2)-1) * size(Img, 1) + interRC(:,1); % see matrix indexing
    
    %tmpSqrtLamda2 = lamda2.^lamda2./abs(lamda1);
    tmpSqrtLamda2 = abs(lamda2).*abs(lamda2)./abs(lamda1);    
    zdot(index) = tmpSqrtLamda2(index);
    %[row, col] = find(lamda2 >= 0);
    %index = (col-1) * size(Img, 1) + row; % see matrix indexing
    %zdot(index) = 0;
    
    %% compute line
    zline = zeros(size(Img, 1), size(Img, 2));
    
    if(brightFlag == 1)
        [row, col] = find(lamda1<0);
    else
        [row, col] = find(lamda1>0);
    end
    
    index = (col-1) * size(Img, 1) + row; % see matrix indexing
    tmpDifLamda = abs(lamda1)-abs(lamda2);
    zline(index) = tmpDifLamda(index);
    
    %% multiplied by sigma*sigma
    zdot = zdot*sigma*sigma;
    zline = zline*sigma*sigma;
   
 if ( i==N_m)
 
       figure,% imshow(mat2gray(zdot))
   %imagesc(mat2gray(zdot));colormap(gray)
   imagesc(zdot);colormap(gray);%title('Enhanced DOT Image');
  s=sprintf('print -depsc %s/zdot_%d_Nm%d,print -djpeg %s/zdot_%d_Nm%d;',subfoldername,iprob,i,subfoldername,iprob,i); eval(s)
 imagesc(zdot);colormap(gray);axis off%title('Enhanced DOT Image');
  s=sprintf('print -depsc %s/nzdot_%d_Nm%d,print -djpeg %s/nzdot_%d_Nm%d;',subfoldername,iprob,i,subfoldername,iprob,i); eval(s)

% % figure,% imshow(mat2gray(zdot))
% %    %imagesc(mat2gray(zdot));colormap(gray)
% %    imagesc(zline);colormap(gray);%title('Enhanced Line ');
% %   s=sprintf('print -depsc %s/zline_%d_Nm%d,print -djpeg %s/zline_%d_Nm%d;',subfoldername,iprob,N_m,subfoldername,iprob,N_m); eval(s)

   
    [filteredImage,filBinImg1,zd_new,zd_new1,terminating_pts1]=fil_Img1(zdot,zline,subfoldername,Img,rad,iprob,N_m);
     zFilterDotImgEach(i, :, :) =filteredImage;
    zd_newEach(i, :, :) = zd_new;   
end
    zDotEach(i, :, :) =zdot;
    zLineEach(i, :, :) = zline;
 
    
    %% compare with the diferent levels and save the largset one
    [row, col] = find(zdot>zdotall);
    index = (col-1) * size(Img, 1) + row; % see matrix indexing
    zdotall(index) = zdot(index);
    %isreal(zdotall)
    
    [row, col] = find(zline>zlineall);
    index = (col-1) * size(Img, 1) + row; % see matrix indexing
    zlineall(index) = zline(index);
    
    clear tmpDifLamda tmpSqrtLamda2 smoothedImg tmpLamda
end
s=sprintf('save %s/T%d.mat dmin dmax N_m...',...
             subfoldername,iprob); eval(s) 


rgbDebug;


