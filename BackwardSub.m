% Author: Nicholas Piercy
% Date : 9/13/2021
% Solve Upper Triangular Matrix Using Backward Substitution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
% U     --> Upper Diagonal Matrix
% b     --> RHS Vector

%OUTPUTS
% x     --> Solution Vector

function [x,check] = BackwardSub(U,b)
x = zeros(length(b),1); % Initialize the solution vector

%last entry into x vector
x(end) = b(end)/U(end,end);
for i = (length(b)-1):-1:1 
    x(i) = (b(i) - (U(i,i+1:end)*x(i+1:end)))/U(i,i);
end
check = U*x;
end