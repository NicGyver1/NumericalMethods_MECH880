% Author: Nicholas Piercy
% Date : 10/18/2021
% Bisection Algorithm for finding a root of a continuous function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% interval  -> Interval which contains a root, [a,b]
% f         -> Function
% tol       -> Error tolerance for convergence
function [c,error,iterationCounter] = bisection(interval,f,tol)
a = interval(1);
b = interval(2);
error = 1; % Initial value of error (dummy)
iterationCounter = 1;
eval_previous = 0; % for Problem 8_HW4 numerical methods UNL 10/19/2021
while error>tol
    c = (a+b)/2; % Midpoint of variable interval [a,b], which contains a root
    message = strcat(num2str(c)," is a root!");
    if f(a)*f(c)<0
        b = c; 
    elseif abs(f(c))<eps
        disp(message)
    else; a = c; 
    end
    iterationCounter = iterationCounter + 1;
    error = (1/(2^(iterationCounter+1)))*(interval(2)-interval(1)); % (1/2^(k+1))*(b-a)
    %Print out iterate information
    eval = double(f(c));
    error_P8_HW4 = abs(eval - eval_previous);
    error = error_P8_HW4;
    if error_P8_HW4<tol
        break
    end
    eval_previous = eval;
    message1 = strcat('x = ',num2str(c),' f(x) = ',num2str(eval),...
        '    |x(k+1)-x(k)| = ',num2str(error));
    disp(message1)
end
disp(message)
message2 = strcat('Iterations to converge = ',num2str(iterationCounter));
message3 = strcat('Absolute error at termination = ',num2str(error));
disp(message2)
disp(message3)
end % End of Bisection Method

