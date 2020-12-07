%% Set parameters and get the needed functions.
D=4.0*10^(-7);
R=0.05;
% dr=0.001;
dr=R/100;
T_0=980;
N=R/dr; % dimension of the mesh. No +1 since dirichlet at r=R.

alpha=3*10^(-4);
% alpha=0.00019845
[A,getbfunc]=getLinearSystem(R, dr, D);

dt=0.1;
T_target=40; %Finished when max temperature of ball is this. Don't change.
t0=0;

%% Get a good estimate of the best cooling function
min_tmax=100; % *Probably* won't find a tmax lower than this.
M=min_tmax/dt;
u_store=zeros(N+1,M); % solution store
tspace=zeros(M,1);

% Initialize u to T_0 everywhere
u_init=zeros(N,1) + T_0;
u = u_init;
u_store(:,1)=[u_init; 980];

U_0=zeros(M,1); % cooling function
U_0(1)=980;

m=2;

du_dr_store=zeros(N,M);
du_dr_store(1)=0;

alpha2=3*10^(-3); % just set to two at the moment, play with this.
cooling_func2=@(t) T_0*exp(-alpha2*t); % change this when you want a new function.

while max(u) > T_target
    guess_max = min(980, U_0(m-1) + 20);
    guess_min = max(20, U_0(m-1) - 100);
    
    if m<30 % let system stabilise slightly before using bisection method.
        bfunc=getbfunc(cooling_func2);
        % pass current time so CrankNicolson knows what to call bfunc with.
        StepCrankNicolson=getCrankNicolson(dt, A, bfunc, tspace(m-1));
        [time, u_next]=StepCrankNicolson(u);
        U_0(m)=cooling_func2(time);
        tspace(m)=time;
        du_dr=d_dr(u_next, dr);
        max_du_dr = max(abs(du_dr));
    else % use bisection method
        while abs(guess_max-guess_min)/2 > 0.0001
            U_0(m)=(guess_max + guess_min)/2;
            if U_0(m-1)<20.1 && U_0(m-1)>20
                U(m)=20;
            end
            % A bit hacky to make the API work.
            % get the corresponding element in vector U_0 for time t.
            get_U_0_value=@(t) U_0(uint64(round((t+dt)/dt)));
            bfunc=getbfunc(get_U_0_value);

            % pass the new bfunc and the current time so
            % that CrankNicolson knows what to call bfunc with.
            StepCrankNicolson=getCrankNicolson(dt, A, bfunc, tspace(m-1));
            [time, u_next]=StepCrankNicolson(u);

            % get gradient via backward difference.
            du_dr=d_dr(u_next, dr);
            max_du_dr = max(abs(du_dr));

            % bisection method
            if max_du_dr > 6000
                guess_min=U_0(m);
            else
                guess_max=U_0(m);
            end
            tspace(m)=time;
        end
    end
    du_dr_store(:,m)=du_dr;
    % store solution and add the boundary condition U_0(t) at r=R
    u_store(:,m) = [u; U_0(m)];
    m=m+1;
    if mod(m,5000)==0 % print progress 
        max(u)
        max_du_dr
    end
    u=u_next;
end

%% Solve again with optimal U_0(t) (averaged to avoid jumps) that was found.
avg_U0 = (U_0(1:end-1)+U_0(2:end))/2;
avg_U0 = [980; avg_U0; 20];
cooling_func_opt=@(t) avg_U0(uint64(round((t+dt)/dt)));

bfunc=getbfunc(cooling_func_opt); % the b(t) term in the linear system Au+b.
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
boundary = getBoundaryFunc(cooling_func_opt);
du_dr_store=zeros(N,1);

while max(u) > T_target
    [time, u] = StepCrankNicolson(u);
    du_dr_store(:,m)=d_dr(u,dr);

    % store solution and add the boundary condition U_0(t) at r=R
    u_store(:,m) = [u; boundary(time)];
    m=m+1;
    tspace(m)=time;
    if mod(m,5000)==0 % print progress 
        max(u)
    end
end

% Add boundary and initial conditions to solutions.
tspace=tspace(1:end-1); % remove last column since u wasnt calculated for that.
rspace=linspace(0,R, R/dr+1);
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
% for my calculated U_0:
%slope around -0.1462 to -0.1465
%T*   = tspace(end) = 7.5419e+03
%T_c* = 6.1331e+03
avg_U0 = (U_0(1:end-1)+U_0(2:end))/2;   

alpha3=1.984498120699633*10^(-4); % just set to two at the moment, play with this.
cooling_func3=@(t) T_0*exp(-alpha3*t); % change this when you want a new function.

figure
hold on
plot(tspace(1:end-1), avg_U0)
% plot(tspace,U_0)
plot(tspace, arrayfun(getBoundaryFunc(cooling_func),tspace))
% plot(tspace, arrayfun(getBoundaryFunc(cooling_func2),tspace))
plot(tspace, arrayfun(getBoundaryFunc(cooling_func3),tspace))

title("Optimal cooling function U_0(t) over time")
lgd=legend("Calculated 'optimal'", "Original exponential, \alpha=3\cdot10^{-4}","Optimal exponential, \alpha=1.98\cdot10^{-4}");
lgd.FontSize = 10;
xlabel("Time")
ylabel("Temperature")

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

