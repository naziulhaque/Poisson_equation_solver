%%Simulation region
% (a,b) -> range of x values
% (c,d) -> range of y values
% (e,f) -> range of z values
% N -> Number of points in any direction   % h -> grid spacing
clc;
clear all;
close all;
a=-5   ; b=5   ;
c=-5   ; d=5   ;
e=-5; f=5;
N=20 ;
h=(b-a)/N;
eps=8.854E-12;

x = linspace(a,b,N);   % grid points x including boundaries
y = linspace(c,d,N);   % grid points y including boundaries
z=linspace(e,f,N);     % grid points z including boundaries
[X,Y,Z] = meshgrid(x,y,z);      % 2d arrays of x,y values

Y=-Y;
Z=-Z;

%%charge distribution (given as input)
    
rho=zeros(N,N,N); 
f=@(x,y,z) (1/(((x)^2+(y)^2+(z)^2)));   % INPUT
for i=1:N
    for j=1:N
        for k=1:N
        rho(i,j,k)=f(X(i,j,k),Y(i,j,k),Z(i,j,k));
        end
    end
end
    
%%Boundary Conditions (dirichlet condition) ->|v <-|^
% V -> Electric Potential 
v=zeros(N,N,N);
 %boundary conditions
v(1, 1:N) = 0;  
v(1:N, N)=  0;
v(N, 1:N)=  0;
v(1:N, 1)=  0;

%%FDM to solve for potential
for iter=1:20
    vi=rho.*(1/6)*(1/eps)*h*h;
for i=2:N-1
    for j=2:N-1
        for k=2:N-1
         vi(i,j,k)=vi(i,j,k)+(1/6)*(v(i-1,j,k)+v(i+1,j,k)+v(i,j-1,k)+v(i,j+1,k)+v(i,j,k-1)+v(i,j,k+1));
        end
    end
end
v=vi;
end

%%Determining the Electric field intensity

ex=zeros(N,N,N);
ey=zeros(N,N,N);
ez=zeros(N,N,N);

for i=2:N-1
    for j=2:N-1
        for k=2:N-1
        ey(i,j,k)=-(v(i,j,k)-v(i+1,j,k))/h;
        ex(i,j,k)=(v(i,j,k)-v(i,j+1,k))/h;
        ez(i,j,k)=-(v(i,j,k)-v(i,j,k+1))/h;
        end
    end
end

for i=2:N-1
    for j=2:N-1
        for k=2:N-1
        ey(i,j,k)=(ey(i,j+1,k)+ey(i,j,k))/2;
        ex(i,j,k)=(ex(i+1,j,k)+ex(i,j,k))/2;
        ez(i,j,k)=(ez(i,j,k+1)+ez(i,j,k))/2;
        end
    end
end


%Normalizing the Electric field
e=zeros(N,N,N);
for i=2:N-1
    for j=2:N-1
        for k=2:N-1
        e(i,j,k)=(ex(i,j,k))^2 + (ey(i,j,k))^2+(ez(i,j,k))^2;
        e(i,j,k)=sqrt(e(i,j,k));
        end
    end
end

for i=2:N-1
    for j=2:N-1
        for j=2:N-1
        if (e(i,j,k)~=0)
            ex(i,j,k)=ex(i,j,k)/e(i,j,k);
            ey(i,j,k)=ey(i,j,k)/e(i,j,k);
            ez(i,j,k)=ez(i,j,k)/e(i,j,k);
        end
        end
    end
end

%%Graphical representation
quiver3(X,Y,Z,ex,ey,ez,'autoscalefactor',3);

