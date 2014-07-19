% Solve wave equation
% using the finite difference method
% written by S J Wiggin (Queen Mary University)
% PDE form: u_{tt}-c^2 u_{xx}=0; c=velocity;

%n=number of points, Nstep=time-step
n = 100; Nstep=100;

% Parameters: L=length of domain;
L=1.0;
velocity=1.0;

%dx and dt<dx/c;
dx = L/(n-1);
dt = 0.9/velocity*dx;
x = 0:dx:L;

% Initial wave profile
q=1;
u=q*exp(-((20/L)*(x-L/2)).^2)';
title_str=...
    'Wave equation: u_{tt}-c^2 u_{xx}=0';

% Store profile for computing and display
v=u;
U=u;
v_initial=u;

%intialization u and v
plot(x,v_initial,'-.','linewidth',2);
axis([0 L -q q]);
set(gca, 'fontsize', 14,'fontweight','bold');
title(title_str);
axis tight;
set(gca,'nextplot','replacechildren');

% start time-stepping

for time=1:Nstep,
    t=time*dt;
    
    % u_{xx}: 2nd derivative at both end
    
    uxx(1)=(v(2)-2*v(1)+v(n))/(dx*dx);
    uxx(2)=(v(3)-2*v(2)+v(1))/(dx*dx);
    uxx(n)=(v(1)-2*v(n)+v(n-1))/(dx*dx);
    uxx(n-1)=(v(n)-2*v(n-1)+v(n-2))/(dx*dx);
    for i=3:n-2,
        uxx(i)=(v(i+1)-2*v(i)+v(i-1))/(dx*dx);
        
    end
     
    %Wave eqution
    % u_{tt}=u_{xx};
    % n^{n+1}=-u^{n-1}+[2u^n+dt*dt*u_{xx}]
    
    u=-u+(2*v+dt*dt*velocity^2*uxx');
    
    % reflecting boundary,
    % which chanmge the p[hase by 180 degrees
    u(1)=0;
    
    %updating arrays by storing results as v
    U=v;
    v=u;
    
    % Plotting the results
    plot(x,v_initial,'-.',x,u,'linewidth',2);
    axis([0 L -q q]);
    picture(time)=getframe;
    
end

%Diaplay movie
disp('playing wave animation ...')
movie(picture,1);