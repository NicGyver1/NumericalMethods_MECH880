% Author: Nicholas Piercy
% Date : 11/1/2021
% Algorithm to add point and compute missing divided difference table
% entries for newtons interpolating polynomial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
% x     --> vector of x values for data points
% y     --> vector of y values for data points

%OUTPUTS
% gamma_new     --> New Coefficient w.r.t. the added point for the Newton Polynomial
% DivDiff_new   --> Expanded Divided Difference Table, of new data point (in an array)
% x_new         --> Total x data points
% y_new         --> Total y data points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gamma_new,DivDiff_new,x_new,y_new] = Newton_DivDiff_PointAdded(x,y,xOriginal_values,DivDiff,gamma)
% Input is added point's (x,y) coordinate
x_new = [xOriginal_values,x];
n_new = length(x_new);
DivDiff_new = DivDiff;
DivDiff_new(n_new,1) = y;
for i = 1:n_new-1 % Columns  
        DivDiff_new(n_new,i+1) = (DivDiff_new(n_new,i)-DivDiff_new(n_new-1,i))/(x_new(n_new)-x_new(n_new-i));
end
% Extract Gamma Values 
% (Coefficients for Newtons Interpolating Polynomial Basis)
gamma_new = vertcat(gamma,DivDiff_new(n_new,n_new));
end  % end of Newton Divided Differences for an Added Point Algorithm