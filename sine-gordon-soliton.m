% Solving the Sine-Gorden equation for soliton propagation
% using the finite differen ce method
%written by S J Wiggin
% PDE form: u_{tt}- u{xx}\alpha*sin(\omega u)

%n=number of points, Nsteps=time-step
n=100; Nstep=100;

%parameters of sine-gordon equation
omega; alpha=15; q=0.1;

% dx and dt<dx;
L=1; dx=L/(n-1);
dt=0.9*dx; x=0:dx:L;

%Initial wave profile