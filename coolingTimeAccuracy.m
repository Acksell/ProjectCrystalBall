clear all;
close all;
R=0.05;
initial_dr = R/10;
initial_dt = 1;

ks = [];
errs = [];

for k = [1 2 3 4 5]
    dr = initial_dr * 2 / (2^k);
    dt = initial_dr * 2 / (2^k);
    [finalTime, finalTemp, u_store, tspace] = getCoolingTime(dt, dr);
    if k > 1
        ks(end+1) = k;
        errs(end+1) = abs(prevFinalTime - finalTime);
    end
    prevFinalTime = finalTime;
end

loglog(ks, errs)

%[t1, u1] = getCoolingTime(dt, dr)
%[t2, u2] = getCoolingTime(dt/2, dr/2)
%diff = (t1-t2)
%diffp = abs(diff)/t1
%%
% 
% maxu = max(u_store, [], 1);
% maxmaxu = max(maxu);
% minmaxu = min(maxu);
% figure
% plot(tspace, -maxu/maxmaxu + minmaxu/maxmaxu + 1)
% figure
% plot(tspace, maxu);
% 
