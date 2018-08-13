function [ph,u]=relax_coarsest(z,ph,beta_h,mu_h,tmarchit,epsilon,dt,lambda1,lambda2,meth)

G=fspecial('gaussian',5,2); % Caussian kernel
Img_smooth=conv2(z,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;nu=100;
g=1./(1+nu*f);  % edge indicator function.
  H=(1+2/pi*atan(ph./epsilon))/2;
%  c1=sum(sum(z.*H))/sum(sum(H));
[vx, vy]=gradient(g);
 A1=sum(sum(H));A2=sum(sum(1-H));nu=0.0001;
for I=1:tmarchit, 
% figure (I)
% mesh(ph)
   H=(1+2/pi*atan(ph./epsilon))/2;
c2=sum(sum(z.*(1-H)))/sum(sum(1-H));  
 c1=sum(sum(z.*H))/sum(sum(H));
Vol=-nu.*( (sum(sum(H))-A1)-(sum(sum(1-H))-A2)  );
     f=lambda1*(z-c1).^2-lambda2*(z-c2).^2; 
    % delta=epsilon./(pi*(epsilon^2+ph.^2));
  f=f+Vol;
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
end

end 
u=c1.*H+c2.*(1-H);

function f = div(nx,ny)
[nxx,junk]=gradient(nx);  
[junk,nyy]=gradient(ny);
f=nxx+nyy;


