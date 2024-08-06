% Author: Nicholas Piercy
% Date : 9/13/2021
% Numerical Methods Homework 2 - Backward Substitution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

A1 = [3 3 1;
     0 -4 -3;
     0 0 2]; % define coefficient matrix
b1 = [12, -10, 4]'; % define loading/RHS vector
[x1,check1] = BackwardSub(A1,b1);

A2 = [3 0 0;
      2 -1 0;
      5 1 -2]; % define coefficient matrix
b2 = [15 10 5]'; % define loading/RHS vector
[x2,check2] = ForwardSub(A2,b2);


% function [x,check] = BackwardSub(A,b)
% 
% x = zeros(length(b),1); % Initialize the solution vector
% 
% %last entry into x vector
% x(end) = b(end)/A(end,end);
% for i = (length(b)-1):-1:1 
%     x(i) = (b(i) - (A(i,i+1:end)*x(i+1:end)))/A(i,i);
% end
% check = A*x;
% end

% function [x,check] = ForwardSub(A,b)
% 
% x = zeros(length(b),1); % Initialize the solution vector
% 
% %first entry into x vector
% x(1) = b(1)/A(1,1);
% 
% for i = 2:length(b) 
%     x(i) = (b(i) - (A(i,1:i-1)*x(1:i-1)))/A(i,i);
% end
% check = A*x;
% end
