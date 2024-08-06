% Author: Nicholas Piercy
% Date : 11/18/2021
% Composite Simpsons Rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%INPUTS
% interval    --> Global Interval [a b]
% myfunc      --> function to be integrated (approximately)
% m           --> Number of subintervals
function [I] = Composite_Simpsons(interval,myfunc,m)
a = interval(1);
b = interval(2);

h = (b-a)/m;
n = 2*m+1; % number of points in each subinterval (2 for Simpsons)
xvals = linspace(a,b,n);
if m == 1
    I = SingleSimpsons(interval,myfunc);
else
    fa = myfunc(a);
    fb = myfunc(b);
    for j = 1:m %first sum term
        f1(j) = myfunc(xvals(2*j));
    end
    f1_term = sum(f1);
    for k = 1:m-1 %second sum term   
        f2(k) = myfunc(xvals(2*k+1));
    end
    f2_term = sum(f2);
    I = (h/6)*(fa + 4*f1_term + 2*f2_term + fb);
end