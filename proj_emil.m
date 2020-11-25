%% Cooling of a crystal ball
clc;clear;close all;
global T_start alpha N D R h dt

% Variables
T_start = 980;                      % Inital temperature
T_end = 40;                         % Ending temperature
T_room = 20;                        % Room temperature
alpha = 3e-4;                       % Cooling rate
D = 4e-7;                           % Thermal diffusivity
R = 5e-2;                           % Radius of the crystal ball
N = 100;                            % Number of space intervals
h=R/N;                              % Step in space
dt = .5;                            % Step in time
T_c = log(T_start/T_room)/alpha;    % Time to reach room temperature

[U, itr] = Crank_Nicolson(T_end);   % Input = desired ending temperature

u_plot(U,itr)
u_inf_plot(U)







% Temperature gradient (?)
M = length(U);
du = [];
for i = 2:M
    du(i) = max(abs(U(3:end,i) - U(1:end-2,i)))/(2*h); 
end
test = max(du);
figure()
plot(linspace(0,dt*M,M),(du))
title('$$ u^{\infty}_r (t) $$','Interpreter','latex')
ylabel('u [°C]')
xlabel('t [s]')


du_max_plot(U) 

%% Functions
function U = U_0(t)
global T_start alpha N
U = ones(N,1)*max(T_start*exp(-alpha*t),20);
end

function A = A_matrix()
global N h R D
% Central difference
% First degree
upper_1 = diag(ones(1,N-1)./(0:h:(R-2*h)),1);
lower_1 = -1*diag(ones(1,N-1)./(h:h:(R-h)),-1);
A_1 = 2*D*(upper_1 + lower_1)/(2*h);

% Second degree
upper_2 = diag(ones(1,N-1),1);
center_2 = diag(-2*ones(1,N));
lower_2 = diag(ones(1,N-1),-1);
A_2 = D*h^-2*(upper_2 + center_2 + lower_2);

A = A_1 + A_2;
A(2,1) = 0;
%BC
A(1,1) = -6*D/(h^2); A(1,2) = 6*D/(h^2);
end

function [U, itr] = Crank_Nicolson(T_stop)
global N dt h D
% Crank-Nicolson
% Pre-Calc
A = A_matrix();
U = zeros(N,N);
U(:,1) = U_0(0);
M1 = sparse(inv(eye(N)-.5*dt*A));
M2 = sparse(eye(N)+.5*dt*A);
b = zeros(N,1); b(end) = 1;

% Method
itr = 1;
while max(U(:,itr)) > T_stop    % Stops at specified temperature
    % Version 1
    U(:,itr + 1) = M1*(M2*U(:,itr) + ...
                    .5*dt*(b.*U_0(dt*itr) + b.*U_0(dt*(itr+1)))*(D/h^2));
    itr = itr + 1;
end

end

function u_plot(u,i)
global R N dt
figure()
mesh(linspace(0,(i)*dt,(i)),linspace(0,R,N),u)
zlabel('u [°C]')
ylabel('r [m]')
xlabel('t [s]')
end

function u_inf_plot(U)
global dt
M = length(U);
u = zeros(1,M);
for i = 1:M
    u(i) = max(U(:,i));
end
figure()
plot(linspace(0,dt*M,M),u)
title('$$ u^{\infty} (t) $$','Interpreter','latex')
ylabel('u [°C]')
xlabel('t [s]')
end

function du_max_plot(U)
global dt N h
M = length(U);
du = zeros(N-1,M);
for i = 2:N-1
    du(i,:) = (abs(U(i+1,:) - U(i,:))); 
end
figure()
plot(linspace(0,dt*M,M),max(du)/(2*h))
title('$$ u^{\infty}_r (t) $$','Interpreter','latex')
ylabel('u [°C]')
xlabel('t [s]')


end
