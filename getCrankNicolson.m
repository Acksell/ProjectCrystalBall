% 1D Crank-Nicolson
function [stepfunc]=getCrankNicolson(dt, A, bfunc)
    dim=length(A); % get dimension for initializing identity matrices.
    
    sparse_left = (speye(dim) - A*dt/2);
    sparse_right = speye(dim) + A*dt/2;
    t=0;
    function [time, u_new]=StepCrankNicolson(u)
        % solve Crank-Nicolson sparse linear system with matlab's backslash.
        u_new = sparse_left\(sparse_right*u + dt*(bfunc(t) + bfunc(t+dt))/2);
        t=t+dt; % increment time.
        time=t; % return time
    end
    stepfunc=@StepCrankNicolson; % return step function
end
