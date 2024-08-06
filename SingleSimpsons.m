% Author: Nicholas Piercy
% Date : 11/17/2021
%Simpsons method over single interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%INPUTS
% interval    --> Global Interval [a b]
% myfunc      --> function to be integrated (approximately)
function [I] = SingleSimpsons(interval,myfunc)
a = interval(1);
b = interval(2);
m = 1; % single interval
h = (b-a)/m;

f1 = myfunc(a);
f2 = myfunc((b+a)/2);
f3 = myfunc(b);
I = (h/6)*(f1 + 4*f2 + f3);

