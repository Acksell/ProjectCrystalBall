%% Set parameters and get the needed functions.
D=4.0*10^(-7);
R=0.05;
% dr=0.001;
dr=R/100;
T_0=980;
N=R/dr; % dimension of the mesh. No +1 since dirichlet at r=R.

alpha=3*10^(-4);
% alpha=0.00019845
cooling_func=@(t) T_0*exp(-alpha*t); % change this when you want a new function.

[A,getbfunc]=getLinearSystem(R, dr, D);

dt=0.1;
T_target=40; %Finished when max temperature of ball is this. Don't change.
t0=0;

%% Solve 
bfunc=getbfunc(cooling_func); % the b(t) term in the linear system Au+b.
StepCrankNicolson = getCrankNicolson(dt, A, bfunc, t0);

% Initialize u to T_0 everywhere
u_init=zeros(N,1) + T_0;
u = u_init;

% Allocate memory for store of solution.
min_tmax=100; % *Probably* won't find a tmax lower than this.
M=min_tmax/dt;
u_store=zeros(N+1,M);
tspace=zeros(M,1);

% Step forward until max temperature is T_target.
m=1;
boundary = getBoundaryFunc(cooling_func);
du_dr_store=zeros(N,1);

while max(u) > T_target
    [time, u] = StepCrankNicolson(u);
    % calculate and store derivative.
    du_dr_store(:,m)=d_dr(u, dr);
    % store solution and add the boundary condition U_0(t) at r=R
    u_store(:,m) = [u; boundary(time)];
    m=m+1;
    tspace(m)=time;
    if mod(m,5000)==0 % print progress 
        max(u)
    end
end

% Add boundary and initial conditions to solutions.
% tspace=tspace(1:end-1); % remove last column since u wasnt calculated for that.
rspace=linspace(0,R, R/dr+1);
u_init=[u_init; boundary(0)];
u_store=[u_init u_store];
du_dr_store=[zeros(N,1) du_dr_store];


%% Plots
figure
% plot(tspace, du_dr)
plot(tspace, max(abs(du_dr_store)))
title("u^{\infty}_r: Max gradient over time", 'FontSize',14)
xlabel("Time",'FontSize',14)
ylabel('du/dr','FontSize',14)
%%
figure
% plot(tspace, du_dr)
rspace_index_used = 1:9:100;
plot(tspace, du_dr_store(rspace_index_used,:))
title("Gradient du/dr over time for different radii", 'FontSize',14)
xlabel("Time",'FontSize',14)
ylabel("du/dr",'FontSize',14)

%%
figure
hold on
% t_to_plot=logspace(0, log(length(tspace)),20);
t_to_plot=linspace(1, length(tspace),100);
for t_i=t_to_plot
    plot(rspace(1:end-1), du_dr_store(:,uint64(floor(t_i))))
end
xlabel("r", 'FontSize',14)
ylabel("du/dr", 'FontSize',14)
title("u_r at different times")


%% 
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

%% 3d plot
figure
optimized_slice = 1:5:length(tspace);
s=surf(tspace(optimized_slice), rspace, u_store(:,optimized_slice)); 
s.EdgeColor = 'flat';
xlabel("Time", 'FontSize',14)
ylabel("r", 'FontSize',14)
zlabel("Temperature", 'FontSize',14)
title("Temperature distribution over time", 'FontSize',14)


