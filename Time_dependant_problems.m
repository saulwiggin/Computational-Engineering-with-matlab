n=100;

L=0.1;
speed=1.0;
m=n-1;
time=1;
dt=L/(n*speed); hh = L/m;
N=time/dt;

for i=1:n,
    x(i)=(i-1)*L/m;
end
x(1)=0; x(n)=L;

for i=1:m,
    E(1,i)=i;
    E(2,i)=i+1;
    h(i)=abs(x(E(2,i))-x(E(1,i)));
end

u=zeros(1,n)'; f=zeros(1,n)';
k=[1 -1; -1 1]; c=[1 1/2; 1/2 1]*1/3;
A=zeros(n,n); C=zeros(n,n);

% C DU/dt+KU=g;
for i=1:m
    A(i,i) = A(i,i)+1/h(i);
    C(i,i)=C(i,i)+h(i)*c(1,1);
    A(i,i+1)=A(i,i+1)-1/h(i);
    C(i,i+1) = C(i,i+1)+h(i)*c(1,2);
    A(i+1,i)=A(i+1,i)-1/h(i);
    C(i+1,i)=C(i+1,i)+h(i)*c(2,1);
    A(i+1,i+1)=A(i+1,i+1)+1/h(i);
    C(i+1,i+1)=C(i+1,i+1)+h(i)*c(2,2);
end

A(n,1)=1 ;A(n,n-1)=0; A(1,1)=1; A(1,2)=0;

C=zeros(n,n);
for i=2:n-1,
    C(i,i)=hh;
end

C(1,1)=hh/2;
C(n,n)=hh/2; Cinv=inv(C);
D=2*eye(n,n)-dt*dt*Cinv*A;

u0=exp(-(40*(x-L/2)).^2)+0.5*exp(-(40*(x-L/4)).^2);
v=u0'; U=v;

plot(x,u0);
axis([0 L -1 1]);
set(gca, 'fontsize', 14, 'fontweight', 'bold');
title('travelling wave');
xlabel('x'); ylabel('u(x,t)');
axis tight;
set (gca, 'nextplot', 'replacechildren');

for t=1:N,
    u=D*U-v;
    v=U;
    U=u;
    plot(x,u,x,u0,'-.','linewidth',2);
    axis([0 L -1 1]);
    V(t)=getframe;
end

movie(V)