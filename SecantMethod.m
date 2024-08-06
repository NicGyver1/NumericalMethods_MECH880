% Author: Nicholas Piercy
% Date : 10/19/2021
% Secant method/iterations for finding a root of a continuous function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% guess     -> initial guesses of root (2 guesses required)
% a         -> slope values to try for a function 
% f         -> Function
% tol       -> Error tolerance for convergence
function [root,errorTracker,numIterations,x] = SecantMethod(guess,f,tol)
error = 1; % Initial value of error (dummy)
eval_previous = error;
j = 1;
message_guess = strcat('Initial guesses = ',num2str(guess(1)),', ',num2str(guess(2)));
disp(message_guess)
while error>tol
    j = j+1;
    if j == 2 
        denominator = double(f(guess(2))-f(guess(1)))/(guess(2)-guess(1));
        x(1) = guess(2);
        x(j) = guess(2) - (double(f(guess(2))/denominator));
    else
        denominator = double((f(x(j-1))-f(x(j-2)))/(x(j-1)-x(j-2)));
        x(j) = x(j-1) - (double(f(x(j-1))/denominator));
    end
    eval = double(f(x(j)));
    format long
    error = abs(eval);
    error_P8_HW4 = abs(eval-eval_previous);
    error = error_P8_HW4;
    message1 = strcat('x = ',num2str(x(j)),' f(x) = ',num2str(eval),...
        '    |x(k+1)-x(k)| = ',num2str(error));
    disp(message1)
    if error_P8_HW4<tol
        break
    end
    eval_previous = eval;   
end
root = x(j); errorTracker = error; numIterations = j;
message = strcat(num2str(root)," is a root!"); disp(message)
message4 = ('--------------------------------');
disp(message4)
end % End of Secant Method

