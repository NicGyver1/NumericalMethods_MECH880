% Author: Nicholas Piercy
% Date: 9/7/2021
% Numerical Methods Homework 1, Computer problem 1.7
clear all, close all, clc

% Variables
x = 1;

format long
%Actual derivative for comparison
fprime_actual = (sec(x))^2;
%approximation of derivative
[finiteDifference_error,centeredDifference_error,h] = fprime(x,fprime_actual);

%plot of magnitidue of error vs h (log log plot)
loglog(h,finiteDifference_error,'k',h,centeredDifference_error,'r')
hold on
legend('Finite Difference Error','Centered Difference Error')
xlabel('h')
ylabel('Error')
title('Computer Problem 1.7')
hold off

function [finiteDifference_error,centeredDifference_error,h] = fprime(x,fprime_actual)
    for k = 0:16
        h(k+1) = 10^(-k);
        finiteDifference_approx = abs((tan(x+h(k+1)) - tan(x))/h(k+1));
        centeredDifference_approx = abs((tan(x+h(k+1)) - tan(x-h(k+1)))/(2*h(k+1)));
        finiteDifference_error(k+1) = abs(finiteDifference_approx-fprime_actual);
        centeredDifference_error(k+1) = abs(centeredDifference_approx-fprime_actual);
    end
end