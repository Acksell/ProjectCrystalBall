% Animated u(r,t) plot

function animated_u(u, R, N, h_t)
% Creates a video of u(r,t) u - y axis , r - x-axis 

% R - sphere radius
% N - intervals
% h_t time step
% t - time
% u - solution to crank nicolson
figure()
myVideo = VideoWriter('myVideoFile'); 
myVideo.FrameRate = 10;  
open(myVideo)
for i=1:40:length(u); % Using step size 40 s to speed it up
    plot(linspace(0,R,N), u(:,i)) %plotting first column from u-matrix over the radius
    seconds = i;
    minutes = ceil(i/60); % To get a rounded number in plot title
    titleText = ['Temperature at time ',num2str(minutes), ' min | Used cooling function: $$U_0 = max(980e^{-\alpha t}, 20)$$'];
    
    title(titleText,'Interpreter','latex' );
    % Making nicer looking axes
    ylim([-10, 1000]); xlim([-0 0.05]);
    xticks([0e-2 1e-2 2e-2 3e-2 4e-2 5e-2]);
    xticklabels({'0' '1' '2' '3' '4' '5'});
    xlabel('Radius [cm]');
    ylabel('Temperature [°C]');
    pause(0.001)
    frame = getframe(gcf); 
    writeVideo(myVideo, frame);
end
close(myVideo)
end