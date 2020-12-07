% Test emil
D=4.0*10^(-7);
R=0.05;
dr=0.001;
T_0=980;
N=R/dr; % dimension of the mesh. No +1 since dirichlet at r=R.

alpha=3*10^(-4); % just set to two at the moment, play with this.
cooling_func=@(t) T_0*exp(-alpha*t); % change this when you want a new function.

[A,bfunc]=getLinearSystem(R, dr, D, cooling_func);

dt=0.1;
StepCrankNicolson = getCrankNicolson(dt, A, bfunc);

% initialize u to T_0 everywhere
u_init=zeros(N,1) + T_0;
u = u_init;

T_target=40; %don't change.

% allocate memoery for store of solution. This might be premature optimisation.
% we start by allocating enough memory to account for rough_tmax/dt steps.
% then the loop will allocate more as it needs. This makes it go a little
% bit faster.
% rough_tmax=0.1;
% M=rough_tmax/dt;
u_store=zeros(N+1,1);
tspace=zeros(10,1);
m=1;
boundary = getBoundaryFunc(cooling_func);
while max(u) > T_target
    [time, u] = StepCrankNicolson(u);
    % store solution and add the boundary condition U_0(t) at r=R
    u_store(:,m) = [u; boundary(time)];
    m=m+1;
    tspace(m)=time;
    if mod(m,5000)==0 % print progress 
        max(u)
    end
end

% add boundary and initial conditions
tspace=tspace(1:end-1); % remove last column since u wasnt calculated for that.
rspace=linspace(0,R, R/dr+1);
tspace=[0; tspace];
u_init=[u_init; boundary(0)];
u_store=[u_init u_store];

% plots
figure
hold on
% t_to_plot=logspace(0, log(length(tspace)),20);
t_to_plot=linspace(1, length(tspace),20);
for t_i=t_to_plot
    plot(rspace, u_store(:,uint64(floor(t_i))))
end
xlabel("r", 'FontSize',14)
ylabel("Temperature", 'FontSize',14)
title("Temperature profiles at different times")

figure
s=surf(tspace, rspace, u_store);
s.EdgeColor = 'flat';
xlabel("Time", 'FontSize',14)
ylabel("r", 'FontSize',14)
zlabel("Temperature", 'FontSize',14)
title("Temperature distribution over time", 'FontSize',14)
%%
animated_u(u_store, R, N+1, dt)
