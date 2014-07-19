% Solving linear elastic problems 
% by the finite element method
% matlab program written
% by S J Wiggin
% usgae : elasticity(ni,nj);
% e.g. elasticity(15,4);

function elasticity(ni,nj)
if nargin <2,
    disp('Usage: elasticity(ni,nj)');
    disp('e.g., elasticity(15,4)');
end
if nargin<1, ni=15; nj=4; end

%triangular mesh/elements
% for 2-D plane stress
% force at top end
fx = 10; fy=0.5;
% applied/fixed displacement
u0=0.0; v0=0.0;

height =2; width=0.4; tunit=0.2;

n=ni*nj; m=2*(ni-1)*(nj-1);

dx=width/(nj-1);
dy=height/(ni-1);

subplot(1,3,1);
pcolor(zeros(ni,nj));
set(gca, 'fontsize',14,'fontweight', 'bold');
title('mesh grid');

A = zeros(2*n,2*n); f=zeros(1,2*n)' ;F=f;

Y=10^6; nu=0.25;

pen=1000*Y;

ushow=1; vshow=1;

k=0;
for i=1:ni,
    for j=1:nj,
        k=k+1;
        x(k)=(j-1)*dx; y(k)=(i-1)*dy;
        node(k)=k;
        xp(i,j)=x(k); yp(i,j)=y(k);
    end
end

k=-1; nk=0;
bcsize=2*(nj+ni-2);

for i=1:ni-1,
    for j=1:nj-1,
        k=k+2;
        nk=node((i-1)*nj+j);
        E(1,k)=nk;
        E(2,k)=nk+1;
        E(3,k)=nj+nk;
        E(1,k+1)=nk+1;
        E(2,k+1)=nj+nk+1;
        E(3,k+1)=nj+nk;
    end
end

%boundary conditions
BC(1:nj)=1:nj;
TOPFORCE(1)=[(ni-1)*nj+1];

area=height*width/m; s=1/(2*area);
G=Y*[1 -nu 0; -nu 1 0; 0 0 (1-nu)/2];
for i=1:m,
    ja = node(E(1,i)); x1=x(ja); y1=y(ja);
    jb=node(E(2,i)); x2=x(jb); y2=y(jb);
    jc=node(E(3,1)); x3=x(jc); y3=y(jc);
    
    %Matrix manipulation
    
    a = ...
        [2*(ja-1)+1 2*ja 2*(jb-1)+1 2*jb 2*(jc-1)+1 2*jc];
    aj = [a;a;a;a;a;a]; ai=aj';
    B=s*[y2-y3 0 y3-y1 0 y1-y2 0;
        0 x3-x2 0 x1-x3 0 x2-x1;
        x3-x2 y2-y3 x1-x3 y3-y1 x2-x1 y1-y2];
    
    %stiffness nmatrix
    K=B'*G*B*tunit*area;
    
    for ii=1:6,
        for jj=1:6,
            A(ai(ii,jj),aj(ii,jj))=...
                K(ii,jj)+A(ai(ii,jj),aj(ii,jj));
        end
    end
end

for i=1:size(BC,2),
    ij=2*BC(i);
    A(ij-1,ij-1)=A(ij-1,ij-1)+pen;
    F(ij-1)=f(ij-1);
    f(ij-1)=f(ij)+pen*u0;
    A(ij,ij)=A(ij,ij)+pen;
    F(ij)=f(ij); f(ij)=f(ij)+pen*v0;
    text(i-0.1,1,'\Delta');
end

for i=1:size(TOPFORCE,2),
    ij=2*TOPFORCE(i);
    f(ij-1)=f(ij-1)+fx;
    F(ij) = f(ij)+fy;
    F(ij)=f(ij);
    text(0.5,ni,'\rightarrow',...
        'fontsize', 20, 'color', 'r');
end

u=A\f;
k=-1;
for i=1:ni,
    for j=1:nj,
        k=k+2;
        U(i,j)=u(k);
        V(i,j)=u(k+1);
        Fx(i,j)=F(k);
        Fy(i,j)=F(k+1);
    end
end

subplot(1,3,2);
pcolor(xp+U*ushow,yp+V*vshow,U);
shading interp; colorbar;
set(gca, 'fontsize', 14, 'fontweight', 'bold');
axis off;
title('U Displacement');
subplot(1,3,3);
pcolor(xp+U*ushow,yp+V*vshow,V);
axis off;
shading interp;
set(gca,'fontsize',14, 'fontweight', 'bold');
colorbar;
title('V displacemnt');