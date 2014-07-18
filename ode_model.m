%-- save this part separetly as ode_model.m
% Belousov-Zhabotinskii reactions
% program by S J Wiggin (Queen Mary University)
% e*dx/dt=q*y-x*y+x*(1-x)
% d*dy/dt=-q*y-x*y+2*p*z
% dz/dt=x-z;
% the following notations are used:
% y(1)=x; y(2)=y; y(3)=z;

% Function for ode_solvers

function dy=ode_model(t,y)
%Initialising 
dy =zeros(3,1);
q=0.005;
p=0.3;
e=0.01;
d=0.15;
% Equations
dy(1)=(y(1)+q*y(2)-y(1)*y(2)-y(1)^2)/e;
dy(2)=(-q*y(2)-y(1)*y(2)+2*p*y(3))/d;
dy(3)=(y(1)-y(3));