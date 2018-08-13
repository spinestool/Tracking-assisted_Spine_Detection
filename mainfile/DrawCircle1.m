function [ph1,ph2]=DrawCircle1(n,m,x10, y10, r,r2)
% Draw a circle on the current figure using ploylines
%
%  DrawCircle(x, y, r, nseg, S)
%  A simple function for drawing a circle on graph.
%
%  INPUT: (x, y, r, nseg, S)
%  x, y:    Center of the circle
%  r:       Radius of the circle

h=1;k=1;ph=zeros(n,m);ph10=zeros(n,m);ph1=zeros(n,m);ph2=zeros(n,m);
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
      ph10(i,j)=-sqrt((xi-x0)^2+(yj-y0)^2)+r2;
    end
 end
a= find(ph>=0);ph1(a)=1;
    b= find(ph10>=0);ph2(b)=1;    
  
end

%phbw=bwdist(ph1);

% figure,imagesc(ph1),colormap(gray)
% hold on
% contour(ph1,[0 0], 'yo','LineWidth',2);  colormap(gray);hold off;
