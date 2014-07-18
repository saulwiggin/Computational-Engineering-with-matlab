%Equations of Belousov-Zhabotinskii reactions
%e*dx/dt=q*y-x*y+x*(1-x)
d*dy/dt=-q*y-x*y+2*p*z
% dz/dt=x-z;
% Save this part as a main function 
% main function
function bz_periodic(time)

if nargin<1,
    time =100;
end
% values between 0.01 and 1.0 are all OK;
yinit = 0.1;
options=odeset('RelTol', 1e-4,...
    'AbsTol', [1e-4 1e-4 1e-5]);
[t,y]=ode15s('ode_model',...
    [0 time],[1.0 1.3 0.5],options);
plot(t,y(:,1),'-.',t,y(:,2),'-',...
    t,y(:,3),'--','linewidth',2);
set(gca, 'fontsize', 14,'fontweight','bold');
legend('x','y','z');