function vbee(n,timesteps)

if nargin<2, timesteps=10; end
if nargin<1, timesteps=20;n=100;end

range=[0 10 0 10];
alpha=0.2;
beta=0.1
dx=(range(2)-range(1))/100;
dy=(range(4)-range(3))/100;
[x,y]=meshgrid(range(1):dx:range(2),...
    range(3):dy:range(4));
z=fun(x,y);

figure(1);
mesh(x,y,z);
view([60 40]);
set(gca, 'fontsize',16,'fontweight','bold');
[xn,yn]=bees(n,range);

figure(2); 
contour(x,y,z,10);hold on;

for i =1:timesteps,
    
    zn=fun(xn,yn);
    zn_max=max(zn);
    xo=xn(zn==zn_max);
    yo=yn(zn==zn_max);
    
    figure(2);
    plot(xn,yn,'.',xo,yo,'d','markersize',10,... 
        'markerfacecolor', 'r');
        axis(range);
    set(gca,'fontsize',16,'fontweight','bold');
    
    figure(3);
    plot(xn,yn,'.',xo,yo,'d',...
        'markersize',10,'markerfacecolor','g');
    axis(range);
    set(gca,'fontsize',16, 'fontweight','bold');
    
    [xn,yn]=bee_moving(xn,yn,xo,yo,alpha,beta);
    
    drawnow;
    
    str=strcat('Current best estimates: xo=',...
        num2str(xo));
    str=strcat(str,'yo = ');
    str=strcat(str,num2str(yo));
    disp(str);
    
end

% initialise zee bees
function [xn,yn]=bees(n,range)
xrange=range(2)-range(1);
yrange=range(4)-range(3);
xn=rand(1,n)*xrange+range(1);
yn=rand(1,n)*yrange+range(3);

function [xn, yn]=bee_moving(xn,yn,xo,yo,alpha,beta)
n=size(yn,2);
xn=xn*(1-beta)+xo*beta+alpha*(rand(1,n)-0.5);
yn=yn*(1-beta)+yo*beta+alpha*(rand(1,n)-0.5);

function [z]=fun(x,y)
z=gauss(x,y,1,1,6)+2*gauss(x,y,2,6,6)...
    +gauss(x,y,1,2,2)+gauss(x,y,1,6,1);
function [z]=gauss(x,y,sigma,mu,nu)
z=exp(-((x-mu).^2+(y-nu).^2)/sigma);
