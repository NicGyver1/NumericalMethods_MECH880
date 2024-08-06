% Author: Nicholas Piercy
% Date : 9/13/2021
% Solve Lower Triangular Matrix Using Forward Substitution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
% L     --> Lower Diagonal Matrix
% b     --> RHS Vector

%OUTPUTS
% x     --> Solution Vector

function [x,check] = ForwardSub(L,b)
x = zeros(length(b),1); % Initialize the solution vector

%first entry into x vector
x(1) = b(1)/L(1,1);
for i = 2:length(b) 
    x(i) = (b(i) - (L(i,1:i-1)*x(1:i-1)))/L(i,i);
end
check = L*x;
end