% Solving wave equation via cellular automata
% written by X S Yang (Cambridge University)
% PDE form: u_{tt}-c^2 u_{xx}=0; c=velocity=1;

clear;
%n=number of points, Nstep=time-step
n =100; Nstep=70;
%Parameters: L=length of domain;
L=n;
% dx=1 and dt=1;
dx = 1; dt=1;
x=1:dx:L;

% Initial wave profile in 2^10=1024 states
q=2^10;
u=fix(q*exp(-((20/L)*(x-L/2)).^2)');
title_str='Cellular automaton: u_{tt}-u_{xx}=0';
% store profile for computing and siaply
v=u; U=u;
v_initial=u;

%initalisation u and v
figure(1);
plot(x,v_initial,'-.', 'linewidth',2);
axis([0 L -q q]);
set(gca, 'fontsize',14,'fontweight','bold');
title(title_str);
axis tight;
set(gca,'nextplot', 'replacechildren');

%Index for spatial loop
I=2:n-1;

% Time-stepping and loop
for time=1:Nstep,
    %periodic at both ends
    uxx(1) =v(2)-2*v(1)+v(n);
    uxx(n)=v(1)-2*v(n)+v(n-1);
    %second derivative
    uxx(I)=v(I+1)-2*v(I)+v(I-1);
    
    %Rules based on: u_tt=u_{xx};
    u = -U++(2*v+uxx');
    u(1)=0; u(n)=0;
    
    U=v;
    v=u;
    plot(x,v_initial,'-.',x,u,'linewidth',2);
   axis([0 L -q q]);
   Ufinal(time, :)=u';
   pause(0.05);
   drawnow;
end

figure(2);
disp('Displaying the ceullular automaton ...');
waterfall(Ufinal);
view([30 60]);
xlabel('u(x,t)');
ylabel('Time (t)');
set(gca,'fontsize',16,'fontweight','bold');
    