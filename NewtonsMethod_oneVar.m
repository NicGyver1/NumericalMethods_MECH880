% Author: Nicholas Piercy
% Date : 10/18/2021
% Newtons method/iterations for finding a root of a continuous function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% guess     -> initial guess of root
% a         -> slope values to try for a function 
% f         -> Function
% tol       -> Error tolerance for convergence

function [root,errorTracker,numIterations,x] = NewtonsMethod_oneVar(guess,f,tol)
fprime = diff(f);
error = 1; % Initial value of error (dummy)
eval_previous = error;
j = 0;
message_guess = strcat('Initial guess = ',num2str(guess));
disp(message_guess)
while error>tol
    j = j+1;
    if j == 1 
        x(j) = guess - (double(f(guess)/fprime(guess)));
    else
        x(j) = x(j-1) - (double(f(x(j-1))/fprime(x(j-1))));
    end
    eval = double(f(x(j)));
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
root = x(j);
errorTracker = error;
numIterations = j;
message = strcat(num2str(root)," is a root!");
disp(message)
message4 = ('--------------------------------');
disp(message4)
end % End of Bisection Method

