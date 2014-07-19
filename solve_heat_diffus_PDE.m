% Solving heat conduction equation
% using the finite differnce method
% written by S J Wiggin (Queen Mary University)
% PDE of the form: u_t=K u_{xx}

%n=number of points, Nsttep=time-step
n=200;
Nstep=200;

%Parameters: L=length of domain
% kappa=heat diffusivity
L=2;
kappa=1.0;

%dx and dt<0.5(dx)^2/kappa
dx=L/(n-1);
dt=(0.46*dx^2)/kappa;
x=0:dx:L;

% Initial temperature profile (square pulse)
q=1;
u=q*(1-abs(sign(x-L*5/8)+sign(x-L*3/8))./2)';

% store profile for computing and display
v=u;
v_initial=u;
title_str='Heat equation: u_t-\kappa u_{xx}=0';

%initalization of u and v
plot(x,v_initial,'-.','linewidth',2);
axis([0 L -q q]);
set(gca, 'fontsize', 14, 'fontweight', 'bold');
title(title_str);
axis tight;
set(gca,'nextplot','replacechildren');

% Time stepping

for time=1:Nstep,
    t=time*dt;
    % u_{xx}: second derivative at bpoth ends
    uxx(1) = (v(2)-2*v(1)+v(n))/(dx*dx);
    uxx(2)=(v(3)-2*v(2)+v(1))/(dx*dx);
    uxx(n)=(v(1)-2*v(n)+v(n-1))/(dx*dx);
    uxx(n-1)=(v(n)-2*v(n-1)+v(n-2))/(dx*dx);
    
    for i=3:n-2,
        uxx(i)=(v(i+1)-2*v(i)+v(i-1))/(dx*dx);
    end
    
    %Heat equation
    u=v+dt*kappa*uxx';
    
    %updating arrays by storing results as v
    v=u;
    
    %updating arrays by storing result as v
    v=u;
    
    %plottinmg the results
    
    plot(x,v_inital,'-.',x,u,'linewidth',2);
    axis([0 L -q q]);
    picture(time)=getframe;
    
end

%display movie

disp('Playing animation ...');
movie(picture,1);