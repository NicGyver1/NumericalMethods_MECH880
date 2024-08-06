% Author: Nicholas Piercy
% Date : 11/18/2021
% Composite Trapezoidal Rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%INPUTS
% interval    --> Global Interval [a b]
% myfunc      --> function to be integrated (approximately)
% m           --> Number of subintervals
function [I] = Composite_Trapezoidal(interval,myfunc,m)
a = interval(1);
b = interval(2);

h = (b-a)/m;
n = m+1; % number of points in each subinterval (2 for mid-pt)
xvals = linspace(a,b,n);

fa = myfunc(a)/2;
fb = myfunc(b)/2;
if m == 1
    I = h*(fa + fb);
else
    for i = 1:n-2
        f(i) = myfunc(xvals(i+1));
    end
    I = h*(fa + sum(f) + fb);
end



