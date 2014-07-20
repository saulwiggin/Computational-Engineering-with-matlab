clear;
n=60;
Nstep=100;

L=n;
kappa=0.25;

dx=1;
x=1:dx:L;

q=2^8
u=fix(q*(1-abs(sign(x-L*5/8)+sign(x-L*3/8))./2)');
v=u;
v_initial=u;
title_str='Heat equation: u_t-\kappa u_{xx}=0';

figure(1);
plot(x,v_initial,'-.','linewidth',2);
axis([0 L -q q]);
set(gca,'fontsize',14,'fontweight','bold');
title(title_str);
axis tight;
set(gca,'nextplot','replacechildren');

I=2:n-1;

for time=1:Nstep,
    
    uxx(1)=v(2)-2*v(1)+v(n);
    uxx(n)=v(1)-2*v(n)+v(n-1);
    uxx(I)=v(I+1)-2*v(I)+v(I-1);
    
    u=v+kappa*uxx';
    
    v=u;
    
    plot(x,v_initial,'-.',x,u,'linewidth',2);
    axis([0 L -q q]);
    Ufinal(time,:)=u';
    pause(0.1);
    drawnow
end

figure(2);
disp('Display results ...');
waterfall(Ufinal);
xlabel('u(x,t)');
ylabel('Time (t)');
set(gca, 'fontsize',16,'fontweight','bold');
view([30 60]);