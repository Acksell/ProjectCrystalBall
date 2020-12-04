clear all;
close all;
R=0.05;
% Largest step sizes to test
initial_dr = R/20;
initial_dt = 1;

% Store dr in each step
drs = [];
% Store dt in each step
dts = [];
% Store power of 2 for each step
ks = [];
% Store difference from last step, for each step
errs = [];

% Halv the step sizes in each iteration and store the difference
for k = [1 2 3 4 5 6 7 8]
    dr = initial_dr * 2 / (2^k);
    dt = initial_dt * 2 / (2^k);
    [finalTime, finalTemp, u_store, tspace] = getCoolingTime(dt, dr);
    if k > 1
        ks(end+1) = k;
        dts(end+1) = dt;
        drs(end+1) = dr;
        errs(end+1) = abs(prevFinalTime - finalTime);
    end
    prevFinalTime = finalTime;
    disp(['Done for: ' num2str(k)])
end

%%

loglog(2./(2.^ks), errs, 'DisplayName', 'Error')
hold on
loglog(2./(2.^ks), (2./(2.^ks)).^2, 'DisplayName', 'y = x^2')

title('Loglog plot of error vs. step size ratio');
xlabel('$\Delta x_k /  \Delta x_{0}$ (and/or) $\Delta t_k / \Delta t_{0}$', 'Interpreter','latex');
ylabel('abs($T_k$ - $T_{k-1})$', 'Interpreter','latex');
legend


%%

plot(dts, errs);
title('Plot of error vs. time step');
xlabel('$\Delta t$', 'Interpreter','latex');
ylabel('Error');

%%

dt = 0.1;
dr = R/100;

[finalTime1] = getCoolingTime(dt, dr)
[finalTime2] = getCoolingTime(dt*2, dr*2)

% dT = finalTime1 - finalTime2 = 1.4838
% time = 11777.51s +- 1.48s
