%% Cooling of a crystal ball: Finding optimal cooling rate, alpha
% Finding the limiting cooling rate alpha* and cooling time T*
clc;clear;close all;
global T_start N D R h dt T_room T_end 

% Variables
T_start = 980;                      % Inital temperature
T_end = 40;                         % Ending temperature
T_room = 20;                        % Room temperature
D = 4e-7;                           % Thermal diffusivity
N = 100;                            % Number of space intervals
h = R/N;                            % Step in space
dt = 10;                            % Step in time

tic
RR = [5 10 15 20]*1e-2;             % Different radius on the crystal ball
RE = find_alpha_star(RR);
toc

% plots
u_inf(RR,RE(:,1))
u_inf_max(RR,RE(:,1))








%% Functions

function S = find_alpha_star(Radius)
global T_end alpha N h t R T_start T_room
S = zeros(2,3);
fprintf('%s %15s %15s %15s', 'Radius [cm]', 'alpha star', 'T star [s]', 'T_c star [s]')
fprintf('\n')
disp('---------------------------------------------------------')
for k = 1:length(Radius)
% Variables
alpha = 1e-4;                       % Cooling rate, start guess
R = Radius(k);
h=R/N;                              % Step in space
t = 0;

U = Crank_Nicolson(T_end, alpha);   % Input = desired ending temperature

du_check = du_max(U);
% Start guesses
z1 = alpha;
z2 = 1e-5;
z = zeros(1,2);
z(1) = z1;
z(2) = z2;

F = zeros(1,2);
F(1) = max(du_check) - 6000;
alpha = z2;
U = Crank_Nicolson(T_end,alpha);   % Input = desired ending temperature
F(2) = max(du_max(U)) - 6000;
eps = 1e-5; % Margin below 6000°C/m


i = 2;
while abs(F(i)) > eps 
    z(i+1) = secant(F(i),F(i-1),z(i),z(i-1));
    alpha = z(i+1);
    U = Crank_Nicolson(T_end,alpha);   % Input = desired ending temperature
    F(i+1) = max(du_max(U))-6000;
    i = i + 1;
    if i > 25
        break
    end
end
S(k,1) = z(end);                        % alpha*
S(k,2) = t;                             % T*
S(k,3) = log(T_start/T_room)/z(end);    % T*_c

fprintf('%5s %20s %15s %15s', num2str(1e2*Radius(k)), num2str(S(k,1)), num2str(S(k,2)), num2str(S(k,3)))
fprintf('\n')
end
fprintf('\n')
end

function z3 = secant(F2,F1,z2,z1)
z3 = z2 - F2*(z2 - z1)/(F2 - F1);
end

function U = Crank_Nicolson(T_stop, alpha)
global N dt t %h D R 
% Crank-Nicolson
% Pre-Calc
A = A_mat();
U = zeros(N,N);
U(:,1) = U_0(0);
M1 = sparse(inv(eye(N)-.5*dt*A));
t = 0;
% Method
itr = 0;
while max(U(:,itr+1)) > T_stop    % Stops at specified temperature
    U(:,itr + 2) = M1*U(:,itr+1) + .5*dt*M1*(A*U(:,itr+1) + b_vec(t,alpha) + b_vec(t+dt,alpha));
    itr = itr + 1;
    t = t + dt;
end

end

function A = A_mat()
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
A(2,1) = 0; % Trouble with A(2,1) = 1e-18
%BC
A(1,1) = -6*D/(h^2); A(1,2) = 6*D/(h^2);
end

function U = U_0(t)
global T_start alpha N
U = ones(N,1)*max(T_start*exp(-alpha*t),20);
end

function b = b_vec(t,a)
global N T_start D h R
b = zeros(N,1);
b(end) = D*(R+h/2)^2/(R^2*h^2)*max(T_start*exp(-a*t),20);
end

function du_max_plot(U)
global dt N h dt
M = length(U);
du = zeros(N-1,M);
for i = 3:N-1
    du(i,:) = (abs(U(i+1,:) - U(i-1,:))); 
end

max_grad = (du_max(U));          % Maximum gradient
ii = 0;
for i = 1:length(max_grad)
   if max_grad(i) > 6000
       ii = i;
       break
   end
end

if ii ~= 0
    plot(linspace(0,dt*M,M),max(du)/(2*h),'linewidth',2)
    title('$$ u^{\infty}_r (t) $$ with respective $$ \alpha^* $$','Interpreter','latex')
    ylabel('u [°C]')
    xlabel('t [s]')
    hold on
else
    plot(linspace(0,dt*M,M),max(du)/(2*h),'linewidth',2)
    title('$$ u^{\infty}_r (t) $$ with respective $$ \alpha^* $$','Interpreter','latex')
    ylabel('u [°C]')
    xlabel('t [s]')
    hold on
end

end

function du_check = du_max(U)
global h %dt t
M = length(U);
du = zeros(1,M-1);
for i = 2:M 
    du(i) = max(abs(U(3:end,i)-U(1:end-2,i)))/(2*h);
end
% figure()
% plot(0:dt:t,du)
% title('checking $$ u^{\infty}_r (t) $$','Interpreter','latex')
% ylabel('u [°C]')
% xlabel('t [s]')
du_check = (du);
end

function u_inf(RR,alpha_s)
global h N T_end R dt t
figure()
for i = 1:length(RR)
    R = RR(i);
    h = R/N;
    alpha = alpha_s(i);
    U = Crank_Nicolson(T_end,alpha);
    M = length(U);
    u = zeros(1,M);
    for j = 1:M
        u(j) = max(U(:,j));
    end
    plot(0:dt:t,u,'linewidth',2)
    hold on
end
title('$$ u^{\infty} (t) $$ with respective $$ \alpha^* $$','Interpreter','latex')
ylabel('u [°C]')
xlabel('t [s]')
legend(['R = ', num2str(RR(1)),' m'],['R = ', num2str(RR(2)),' m'],...
    ['R = ', num2str(RR(3)),' m'],['R = ', num2str(RR(4)),' m'])
end

function u_inf_max(RR,alpha_s)
global h N T_end R
figure()
for i = 1:length(RR)
    R = RR(i);
    h = R/N;
    alpha = alpha_s(i);
    U = Crank_Nicolson(T_end,alpha);
    M = length(U);
    u = zeros(1,M);
    for j = 1:M
        u(j) = max(U(:,j));
    end
    du_max_plot(U)
end
ylabel('u [°C]')
xlabel('t [s]')
legend(['R = ', num2str(RR(1)),' m'],['R = ', num2str(RR(2)),' m'],...
    ['R = ', num2str(RR(3)),' m'],['R = ', num2str(RR(4)),' m'])
end
