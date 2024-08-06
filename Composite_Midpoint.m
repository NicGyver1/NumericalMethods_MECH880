% Author: Nicholas Piercy
% Date : 11/18/2021
% Composite Midpoint Rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%INPUTS
% interval    --> Global Interval [a b]
% myfunc      --> function to be integrated (approximately)
% m           --> Number of subintervals
function [I] = Composite_Midpoint(interval,myfunc,m)
a = interval(1);
b = interval(2);

h = (b-a)/m;
n = m+1; % number of points in each subinterval (2 for mid-pt)
xvals = linspace(a,b,n);
I = 0; % initialize the soln sum
for i = 1:m
    arg = (xvals(i)+xvals(i+1))/2;
    I = I + myfunc(arg);
end
I = h*I;