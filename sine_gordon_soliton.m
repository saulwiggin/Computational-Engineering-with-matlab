% Solving the Sine-Gorden equation for soliton propagation
% using the finite differen ce method
%written by S J Wiggin
% PDE form: u_{tt}- u{xx}\alpha*sin(\omega u)

%n=number of points, Nsteps=time-step
n=100; Nstep=100;

%parameters of sine-gordon equation
omega=20; alpha=15; q=0.1;

% dx and dt<dx;
L=1; dx=L/(n-1);
dt=0.9*dx; x=0:dx:L;

%Initial wave profile
u=q*exp(-((20/L)*(x-L/2)).^2');
title_str = 'u_{tt}-u_{xx}=\alpha sin(\omega u))';

%iNITIALISATION of arrays for later use
U=u; v=u; v_inital=u;
set(gca,'fontsize',14,'fontweight','bold');
plot(x,v_inital,'-.','linewidth',2);
axis([0 L -q q])
title(title_str);
axis tight;
set(gca, 'nextplot','replacechildren');

for time =1:Nstep,
    t=time*dt;
    % u_{xx}: 2nd derivative/period at both ends
    uxx(1) =(v(2)-2*v(1)+v(n))/(dx*dx);
    uxx(2)=(v(3)-2*v(2)+v(1))/(dx*dx);
    uxx(n)=(v(1)-2*v(n)+v(n-1))/(dx*dx);
    uxx(n-1)=(v(n)-2*v(n-1)+v(n-2))/(dx*dx);
    
    for i=3:n-2,
        uxx(i)=(v(i+1)-2*v(i)+v(i-1))/(dx*dx);
    end
    
    %Sine Gordon Equation:
    % u_{tt}=u_{xx}-slpah sin(omega*u)
    % or n^{n+1}=-u^{n-1}+[2u^n...
    % +dt*dt(u_{xx}-a sin(omega*u))]
    
    u=-U+(2*v+dt*dt*(uxx'-alpha*sin(omega*v)));
    
    % Reflecting boundary by fixing u(1)=u(n)=0
    u(1)=0; u(n)=0;
    
    %updating the arrays
    U=v; v=u;
    %Plotting the results
    plot(x,v_inital,'-.',x,u,'linewidth',2);
    axis([0 L -q q]);
    picture(time)=getframe;
end
% Play movie

disp('Playing the animation ...');
movie(picture,1);
