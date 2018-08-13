function [ph,phbw,ph1]=DrawNewROI_circle(n,m,centroids, r)
x10=centroids(:,1);y10=centroids(:,2);

h=1;k=1;ph=zeros(n,m);ph1=zeros(n,m);
% [xi,yj]=meshgrid(n,m); 
[l,s]=size(x10);
for q=1:l 
   
 y0=x10(q); x0=y10(q);
 
 for i=1:n
    xi=(i-1/2)*h;
    for j=1:m
        yj=(j-1/2)*k;
       ph(i,j)=-sqrt((xi-x0)^2+(yj-y0)^2)+r;
      %  ph(i,j)=abs(xi-x0)+abs(yj-y0)-r1;
     
    end
 end
 
a= find(ph>=0);ph1(a)=1;
   
 
end

phbw=bwdist(ph);

% figure,imagesc(ph1),colormap(gray)
% hold on
% contour(ph1,[0 0], 'yo','LineWidth',2);  colormap(gray);hold off;
