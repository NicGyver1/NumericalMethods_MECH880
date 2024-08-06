% Author: Nicholas Piercy
% Date: 9/7/2021
% Numerical Methods Homework 1 Computer Problem 1.10
clear all, close all, clc
%input Variables as vector
a = [6,6*10^(154),0,1,1,10^(-155)];
b = [5,5*10^(154),1,-10^(5),-4,-10^(155)];
c = [-4,-4*10^(154),1,1,3.999999,10^(155)];
% Function call and print roots of corresponding coefficients
for i = 1:length(a)
    [xint_1(i),xint_2(i)] = quadraticEQ(a(i),b(i),c(i));   
    fprintf('%5.5g	%5.5g	%5.5g	%5.5g	%5.5g\n',[xint_1(i);xint_2(i);a(i);b(i);c(i)])
end
function [xint_1,xint_2] = quadraticEQ(a,b,c)
%minimize overflow by dividing all terms by largest numerator term
denominator = max(abs([a,b,c]));
a = a/denominator; b = b/denominator; c = c/denominator;
complex_test = b^2 - 4*a*c;
if a==0 % unusual input a = 0
    xint_1 = -c/b; xint_2 = NaN;
elseif complex_test < 0 % complex roots
    xint_1 = 'imaginary'; xint_2 = 'imaginary';
elseif (complex_test)^(1/2) ==0 %if cancellation inside sqrt ==> same x-intercepts/roots
    xint_1 = (-b/(2*a));
    xint_2 = (-b/(2*a));
else %try to avoid catastrophic cancellation
    if b>0 % avoids catastrophic cancellation by subtracting sqrt from a negative term(-b)
        xint_1 = (-b-(complex_test)^(1/2))/(2*a); xint_2 = (2*c)/(-b-(complex_test)^(1/2));
    else % avoids catastrophic cancellation by adding sqrt to a positive term(-b)
        xint_1 = (-b+(complex_test)^(1/2))/(2*a); xint_2 = (2*c)/(-b+(complex_test)^(1/2));
    end
end
end % end of quadraticEQ function 
    