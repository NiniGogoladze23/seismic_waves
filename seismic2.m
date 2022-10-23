clear all
clc
%constants
Lx=204;
Lz=204;
T=150;
nx=52;  
nz=52;                                                  
dx=Lx/(nx-1);  
dz=Lz/(nz-1);
dt=0.5; 
nt=round(T/dt);
x=linspace(0,Lx,nx); 
z=linspace(0,Lz,nz);
t=dt.*linspace(0,nt,nt);
%Cz=Cx;
Cx=1;
Cxp=0.0001;
gammax= Cx*dt/dx;
gammaxp= Cxp*dt/dx;
loc=7;
A=2;

g=zeros(nx-1,nz-1);
for i=1:nz
    if i<=(nz-1)/3;
        g(i,:)=gammax;
    elseif (i>(nz-1)/3 & i<=2*(nz-1)/3)
        g(i,:)=gammaxp;
    else
        g(i,:)=gammax;
    end
end
p=zeros(nx,nz,nt);
for i = 1:nx
    for k = 1:nz
        p(i, k, 1) = -A*exp(-((dx*i-Lx/2)/loc)^2-((dz*k-Lz/2)/loc)^2);
    end
end

for n = 2:nt
for i = 2:nx-1
    for k = 2:nz-1
p(i,k,n+1)= 2*p(i,k,n)-p(i,k,n-1)+gammax^2*(p(i+1,k,n)-2*p(i,k,n)+p(i-1,k,n))+g(i,k)^2*(p(i,k+1,n)-2*p(i,k,n)+p(i,k-1,n));%+(dt)^2*100*exp(-n*((i^2+k^2)/2));
    end
end
im=image(x,z,p(:,:,n));
colorbar;
xlabel('x');
ylabel('z');
pause(0.001)
end

