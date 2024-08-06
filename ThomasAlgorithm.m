% Author: Nicholas Piercy
% Date : 9/19/2021
% Thomas Algorithm for Tridiagonal Systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Thom_U,Thom_L] = ThomasAlgorithm(A)
n = length(A(:,1));
Thom_U = zeros(n,n);
Thom_L = eye(n);

Thom_U(1,1) = A(1,1);
for i = 1:n-1
    Thom_U(i,i+1) = A(i,i+1);
end
for i = 2:n
    Thom_L(i,i-1) = A(i,i-1)/Thom_U(i-1,i-1);
    Thom_U(i,i) = A(i,i) - Thom_L(i,i-1)*Thom_U(i-1,i);
end
end

