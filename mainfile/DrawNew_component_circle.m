function [ph,ph1]=DrawNew_component_circle(n,m,centroids, r)
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
     
     
    end
 end
 
a= find(ph>=0);ph1(a)=1;
   
 
end



