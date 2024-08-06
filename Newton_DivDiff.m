% Author: Nicholas Piercy
% Date : 11/1/2021
% Newton Divided Differences Algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
% x     --> vector of x values for data points
% y     --> vector of y values for data points

%OUTPUTS
% gamma     --> Coefficients Of the Newton Polynomial
% DivDiff   --> Divided Difference Table (in an array)
% n         --> number of data points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gamma,DivDiff,n] = Newton_DivDiff(x,y)
n = length(x); % number of data points
gamma = zeros(n,1); % instantiating vector to store gamma values

DivDiff = zeros(n,n);
DivDiff(:,1) = y(:);
for i = 1:n-1 % Columns
    for k = 1:(n-i) % Rows
        DivDiff(k+i,i+1) = (DivDiff(k+i,i)-DivDiff(k+i-1,i))/(x(k+i)-x(k));    
    end
end
% Extract Gamma Values 
% (Coefficients for Newtons Interpolating Polynomial Basis)
for i = 1:n
    gamma (i) = DivDiff(i,i);
end
end % end of Newton Divided Differences Algorithm













