function [coolingTime, finalTemp, u_store, tspace, rspace] = getCoolingTime(dt, dr)

D=4.0*10^(-7);
R=0.05;
T_0=980;
N=R/dr; % dimension of the mesh. No +1 since dirichlet at r=R.
t0=0;

alpha=3*10^(-4); % just set to two at the moment, play with this.
cooling_func=@(t) T_0*exp(-alpha*t); % change this when you want a new function.

[A,getbfunc]=getLinearSystem(R, dr, D);
bfunc=getbfunc(cooling_func);
StepCrankNicolson = getCrankNicolson(dt, A, bfunc, t0);

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
    if m > 1
        prev_time = time;
        prev_max_u = max(u);
    end
    [time, u] = StepCrankNicolson(u);
    % store solution and add the boundary condition U_0(t) at r=R
    % Don't store values unless needed (to run faster)
    %u_store(:,m) = [u; boundary(time)];
    tspace(m)=time;
    m=m+1;
    if mod(m,50000)==0 % print progress 
        max(u)
    end
end

% Interpolate last step linearly to reach the target temperature exactly
last_step_slope = (max(u) - prev_max_u) / dt;
time_diff_to_target = (T_target - prev_max_u)/last_step_slope;

rspace=linspace(0,R, R/dr+1);
tspace=[0; tspace];
u_init=[u_init; boundary(0)];
u_store=[u_init u_store];

% Final cooling time to target temperature
coolingTime = prev_time + time_diff_to_target;
% Temperature of last step (not target temperature)
finalTemp = max(u);

end

