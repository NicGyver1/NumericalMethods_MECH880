% Author: Nicholas Piercy
% Date : 10/18/2021
% Newtons method/iterations for finding a root of a continuous function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% A     -> Overdetermined Matrix A (M>n --> more rows than columns)
% b     -> slope values to try for a function 
% OUTPUTS
% M     -> LHS matrix normal equations
% d     -> Corresponding RHS normal equations

function [M,d] = NormalEqs(A,b)
M = A'*A;
d = A'*b;
end