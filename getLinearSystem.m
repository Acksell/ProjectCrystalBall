% Gets the matrix A and the function b(t) for the linear system
% du/dt = Au + b(t)
% In this case the system is a crystal ball being cooled.
function [A, getbfunc]=getLinearSystem(R, dr, D)
    N=R/dr; % between dr and R-dr, add r=0 to diags after.
    rspace=linspace(0, R-dr, N);

    subdiag = arrayfun(@(r_i) ((r_i-dr/2)/r_i)^2, rspace(2:end));
    maindiag= arrayfun(@(r_i) -( (r_i+dr/2)^2 + (r_i-dr/2)^2 )/((r_i)^2), rspace);
    supdiag = arrayfun(@(r_i) ((r_i+dr/2)/r_i)^2, rspace(2:end-1));
    subdiag=[subdiag 0]; % extend by 1 to match dimension, value doesnt matter.
    supdiag=[0 0 supdiag]; % extend by 2 to match dimension, value doesnt matter.
    % set first row, this implements the neumann condition at r=0.
    maindiag(1) = -6;
    supdiag(2) = 6;
    
    % Generate initial A matrix and then impose neumann boundary condition.
    A=spdiags([subdiag' maindiag' supdiag'],-1:1,N,N);
    %  impose neumann at r=0, see derivations for explanation.
%     A(1,1)=-6;
%     A(1,2)=6;
    A=D*A/(dr^2); % scale by the common factors.
    

    % return a function which abstracts the translation from cooling
    % function to boundary condition to the b term in the linear system.
    function [bfunc]=getbfunction(cooling_func)
        % Dirichlet at r=R is responsible for the b-term in the linear system.
        boundary_function = getBoundaryFunc(cooling_func);
        function [b]=b(t)
            b=zeros(N,1);
            b(end)=D*boundary_function(t)*(R+dr/2)^2/(R^2*dr^2);
        end
        bfunc=@b;
    end
    getbfunc=@getbfunction;
end