function [filteredImage,filBinImg1,zd_new1,zd_new2,terminating_pts1]=fil_Img1(zd,zl,subfoldername,f,rad,iprob,N)
 global ffb zd_treshhold skeleton_tr rad2

filteredImage= im2bw(zd, zd_treshhold);%filteredImage= im2bw(zd, 0.5);
figure
imagesc(filteredImage),%title('Enhanced DOT Image with graythresh');%colormap(gray)
% set(t,'color','r','fontsize',12)
colormap(gray)
s=sprintf('print -depsc %s/zdotgraythresh_%d_Nm%d,print -djpeg %s/zdotgraythresh_%d_Nm%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)

figure;imagesc(filteredImage);colormap(gray);axis off
s=sprintf('print -depsc %s/nzdotgraythresh_%d_Nm%d,print -djpeg %s/nzdotgraythresh_%d_Nm%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)

%g1 = bwmorph(filteredImage, 'thin', 5);
%filteredLine = im2bw(zl, graythresh(zl));

fa=f./max(max(f));
filteredLine = im2bw(fa, skeleton_tr);%%%%%%%%%%%%%
% filteredLine = im2bw(fa, graythresh(fa));%%%%%%%%%%%%%
%filteredLine = im2bw(zl, 0.9);%%%%%%%%%%%%%
g2 = bwmorph(filteredLine, 'thin', Inf);

g2(1:end,1:ffb)=0;g2(1:ffb,1:end)=0;
g2(end-ffb:end,1:end)=0;g2(1:end,end-ffb:end)=0;
% figure,imagesc(g2),colormap(gray)
terminating_pts = find_skel_ends(g2,'testing');

figure;
    imagesc(g2);colormap(gray)
	hold on;
	plot(terminating_pts(:,1),terminating_pts(:,2),'r*');
    %title('Red points indicated detected terminal points');
    s=sprintf('print -depsc %s/skeleton_%d_Nm%d,print -djpeg %s/skeleton_%d_Nm%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)

figure;imagesc(g2);colormap(gray);hold on;
	plot(terminating_pts(:,1),terminating_pts(:,2),'r*');axis off
    s=sprintf('print -depsc %s/nskeleton_%d_Nm%d,print -djpeg %s/nskeleton_%d_Nm%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)


%terminating_pts1=putout_ex( zd,terminating_pts);
% terminating_pts1=putout_ex1( filteredImage,terminating_pts);
terminating_pts1=putout_ex2( filteredImage,terminating_pts);
 %terminating_pts1=putout_ex3( filteredImage,terminating_pts);

figure,imagesc(zd),hold on,plot(terminating_pts1(:,1),terminating_pts1(:,2),'r*');colormap(gray)
 s=sprintf('print -depsc %s/dotendp_%d_Nm%d,print -djpeg %s/dotendp_%d_Nm%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)

 figure,imagesc(zd),hold on,plot(terminating_pts1(:,1),terminating_pts1(:,2),'r*');colormap(gray);axis off
 s=sprintf('print -depsc %s/ndotendp_%d_Nm%d,print -djpeg %s/ndotendp_%d_Nm%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)

 [ph1,ph2]=DrawCircle_on_t(zd,terminating_pts1,rad,rad2);
% zd_new1=zd.*ph1;zd_new2=zd.*ph2;%%%using zdot
zd_new1=f.*ph1;zd_new2=f.*ph2;%%%%using the image
%zd_2=zd.*ph10;%zd_2=double(I).*ph10;
%zd_new1=ph1; Same
%figure,imagesc(zd_new),hold on,plot(terminating_pts1(:,1),terminating_pts1(:,2),'r*');
figure, imagesc(f),hold on; colormap(gray);contour(zd_new1,[0 0], 'ro');
s=sprintf('print -depsc %s/cont_%d_Nm%d,print -djpeg %s/cont_%d_Nm%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)

figure, imagesc(f),hold on; colormap(gray);contour(zd_new1,[0 0], 'ro');axis off
s=sprintf('print -depsc %s/ncont_%d_Nm%d,print -djpeg %s/ncont_%d_Nm%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)

% figure, imagesc(f),hold on; colormap(gray);contour(zd_2,[0 0], 'ro');
% s=sprintf('print -depsc %s/cont_%d2_Nm%d,print -djpeg %s/cont_%d2_Nm%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)


% filteredImage = bwmorph(filteredImage, 'open');
% filteredImage = bwmorph(filteredImage, 'open');
% filteredImage = bwmorph(filteredImage, 'open');

 filteredImage = bwmorph(zd_new1, 'open');
  filBinImg1= bwmorph(zd_new2, 'open');
%filteredImage=zd_new1;
figure, imagesc(filteredImage);colormap(gray)
s=sprintf('print -depsc %s/bwmorph_%d_Nm%d,print -djpeg %s/bwmorph_%d_Nm%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)
% figure, imagesc(filBinImg1);colormap(gray)
% s=sprintf('print -depsc %s/bwmorph_%d_Nm%d,print -djpeg %s/bwmorph_%d_Nm%d;',subfoldername,iprob,N,subfoldername,iprob,N); eval(s)
