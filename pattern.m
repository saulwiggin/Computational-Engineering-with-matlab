% Pattern formation: a 20 line matlab program
% PDE form: u_t=D*(u_{xx}+u_{yy}+gamma*q(u)
% written by S J Wiggin (Queen Mary university)
% Usage: pattern(200)

function pattern(time)
%input number of time steps
if nargin<1, time=100; end

%initalise parametes n=100;
% D=0.2; gamma=0.5;
n=100; D=0.2; gamma=0.5;

%non-linear partial differential equation
% the solution of the PDE is obtained by
% finite difference method
% assuming dx=dt=-dt=1
% u_t=D(u_{xx}+u_{yy})+gamma*q(u);
%qstr='u.*(1-u)';

qstr='u.*(1-u)';

%set inital values of u randomly
u=rand(n,n); grad=u*0;

%index for u(i,j) and loop
I=2:n-1; J=2:n-1;

%Time stepping
for step=1:time,
    %laplace gradient/updating the equation
    grad(I,J) =u(I,J-1)+u(I,J+1)+u(I-1,J)+u(I+1,J);
    u=(1-4*D)*u+D*grad+gamma*eval(qstr);
    % show results
    pcolor(u);
    shading interp;
    % coloring and colorbar
    h = colorbar; colormap jet;
    disp (strcat('Time step = ', num2str(step)));
    drawnow;
    set([gca h], 'fontsize', 14, 'fontweight', 'bold');
end

%plot as a surface
figure(2);
surf(u); shading interp;
set(gca, 'fontsize', 14, 'fontweight', 'bold');
view([-25 70]);