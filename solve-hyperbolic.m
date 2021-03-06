% Solving the 1-D hyperbolic equation
% using the finite difference method
% written by S J Wiggin (Queen Mary University)
% PDE form: u_t + c*u_{x}=0;

%n=number of points, Nstep=time-step
n=100;  Nstep=250;

%Parameters of 1st-order hyperbolic equations
velocity=1.0;
L=1; dx=L/(n-1);
dt=0.1*dx; x=0:dx:L;
% Initial wave profile
q=0.1;
u=q*(0.5*exp(-((20/L)*(x-L/4)).^2)...
    +exp(-((20/L)*(x-L/2)).^2))';
title_str='Hyperbolic equation: u_t+c u_x=0';

%Intialization
U=u; v=u; v_initial=u;

plot(x,v_inital,'-','linewidth',2);
axis([0 L -q q]);
set(gca, 'fontsize',14,'fontweight','bold');
title(title_str);
axis tight;
set(gca, 'nextplot','replacechildren');

for time=1:Nstep,
    t=time*dt;
    
    %u_x: First derivative & periodic boundaries
    ux(1)=(v(2)-v(n))/(2*dx);
    ux(2)=(v(3)-v(1))/(2*dx);
    ux(n)=(v(1)-v(n-1))/(2*dx);
    
    for i=3:n-2,
        ux(i)=(v(i+1)-v(i-1))/(2*dx);
    end
    
    %First-order hyperbolic equation
    % u_t+gamma u_x=0
    % using Staggered Leapfrog method
    
    u = U-velocity*dt*ux;
    
    %Updating the arrays
    U=u; v=u;
    % Plotting the results
    plot(x,v_initial,'-.',x,u,'inewidth',2);
    axis([0 L -q q]); pause(0.01);
    drawnow;
end
