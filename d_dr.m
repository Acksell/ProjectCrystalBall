function [du_dr]=d_dr(u, dr)
    du_dr=zeros(length(u),1); % also sets du_dr(1)=0 because of B.C.
    for j=2:(length(u)) %backward diff. We get 49 instead of 48 points.
        du_dr(j)=(u(j)-u(j-1))/(dr);
    end
end