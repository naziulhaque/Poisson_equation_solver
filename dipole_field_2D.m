%%Simulation region
% (a,b) -> co-ordinate of lower left corner of the square region
% (c,d) -> co-ordinate of upper right corner of the square region
% (c-a)==(d-b)
% N -> Number of points in any direction   % h -> grid spacing
clc;
clear all;
close all;
a=-10   ; b=-10   ;
c=10   ; d=10   ;
N=200 ;
h=(c-a)/N;
eps=8.854E-12;

x = linspace(a,c,N);   % grid points x including boundaries
y = linspace(b,d,N);   % grid points y including boundaries

[X,Y] = meshgrid(x,y);      % 2d arrays of x,y values

Y=-Y;


%%charge distribution
%Point Charge
% Location of point charge is <qx>, <qy> | Value of charge is <q>   (INPUT)
    
rho=zeros(N,N); 
qx=[ -1 1];        % INPUT
qy=[ 0 0 ];        % INPUT
q =[2 -2].*(1E-6); % INPUT

for k=1:length(q)
    qi(k)=(d-qy(k))/h; qi(k)=floor(qi(k))+1;
    qj(k)=(qx(k)-a)/h;  qj(k)= floor(qj(k))+1;
    rho(qi(k),qj(k))=rho(qi(k),qj(k))+q(k)/h/h;       
end
    

%%Boundary Conditions (dirichlet condition) ->|v <-|^
% V -> Electric Potential 
v=zeros(N,N);
v(1, 1:N) = 0;  
v(1:N, N)=  0;
v(N, 1:N)=  0;
v(1:N, 1)=  0;

%%FDM to solve for potential
for iter=1:200
    vi=rho.*(1/4)*(1/eps)*h*h;
for i=2:N-1
    for j=2:N-1
        vi(i,j)=vi(i,j)+(1/4)*(v(i-1,j)+v(i+1,j)+v(i,j-1)+v(i,j+1));
    end
end
v=vi;
end

%%Determining the Electric field intensity

ex=zeros(N,N);
ey=zeros(N,N);

for i=2:N-1
    for j=2:N-1
        ey(i,j)=-(v(i,j)-v(i+1,j))/h;
        ex(i,j)=(v(i,j)-v(i,j+1))/h;
    end
end

for i=2:N-1
    for j=2:N-1
        ey(i,j)=(ey(i,j+1)+ey(i,j))/2;
        ex(i,j)=(ex(i+1,j)+ex(i,j))/2;
    end
end


%Normalizing the Electric field
e=zeros(N,N);
for i=2:N-1
    for j=2:N-1
        e(i,j)=(ex(i,j))^2 + (ey(i,j))^2;
        e(i,j)=sqrt(e(i,j));
    end
end

for i=2:N-1
    for j=2:N-1
        if (e(i,j)~=0)
            ex(i,j)=ex(i,j)/e(i,j);
            ey(i,j)=ey(i,j)/e(i,j);
        end
    end
end

%%Graphical representation

scatter(diag(X(qi,qj))'-h/2,diag(Y(qi,qj))'+h/2,100,'filled');
hold on
quiver(X,Y,ex,ey,'autoscalefactor',0.7)
figure 
surf(X,Y,v);
figure
s=0.7;
for i=a:s:c
streamline(X,Y,ex,ey,a:s:c, i.*ones(1,length(a:s:c)));
hold on
end