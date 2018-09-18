function [ph,u,c1,c2]=relax_coarsest2(z,ph,beta_h,mu_h,tmarchit,epsilon,dt,lambda1,lambda2,meth,nu)

G=fspecial('gaussian',5,2); % Gaussian kernel
Img_smooth=conv2(z,G,'same');  % smooth image by Gaussian convolution
%ph=bwdist(ph); %ph=-double(ph);
 %Img_smooth = z;
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2; %nu=100;
g=1./(1+nu*f);  % edge indicator function.

[vx, vy]=gradient(g);
for I=1:tmarchit, 
% figure (I)
% mesh(ph)
   H=(1+2/pi*atan(ph./epsilon))/2;
  % H=Heaviside1(ph,epsilon); 

c1=sum(sum(z.*H))/sum(sum(H));c2=sum(sum(z.*(1-H)))/sum(sum(1-H));  

     f=lambda1*(z-c1).^2-lambda2*(z-c2).^2; 
    % delta=epsilon./(pi*(epsilon^2+ph.^2));
  
delta=epsilon./(pi*(epsilon^2+ph.^2));

 [phi_x,phi_y]=gradient(ph); dom=sqrt(phi_x.^2 + phi_y.^2);
    Nx=phi_x./(dom+beta_h); % add a small positive number to avoid division by zero
    Ny=phi_y./(dom+beta_h);
    divK=div(Nx,Ny);

if meth==1
ph=ph + dt*delta.*(mu_h.*divK -f);
else
edgeTerm=(vx.*Nx+vy.*Ny) + g.*divK;
%edgeTerm=delta.*(vx.*Nx+vy.*Ny) + delta.*g.*divK;
ph=ph + dt*delta.*(mu_h.*edgeTerm -f);
% ph=bwdist(ph); ph=double(ph);
end

end 

%if c1<c2;ph=-ph;end
u=c1.*H+c2.*(1-H);

function f = div(nx,ny)
[nxx,junk]=gradient(nx);  
[junk,nyy]=gradient(ny);
f=nxx+nyy;


