% Author: Nicholas Piercy
% Date : 10/18/2021
% Newtons method/iterations for finding a root of a continuous function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% guess     -> initial guess of root
% a         -> slope values to try for a function 
% f         -> Function
% tol       -> Error tolerance for convergence

function [root,errorTracker,numIterations,x] = NewtonsMethod(guess,a,f,tol)
fprime = diff(f);
for i = 1:length(a) % to try all values of a
    error = 1; % Initial value of error (dummy)
    j = 0;
    message_guess = strcat('Initial guess = ',num2str(guess(i)));
    disp(message_guess)
    while error>tol
        j = j+1;
        if j == 1 
            x(j,i) = guess(i) - (double(f(guess(i),a(i))/fprime(guess(i),a(i))));
        else
            x(j,i) = x(j-1,i) - (double(f(x(j-1,i),a(i))/fprime(x(j-1,i),a(i))));
        end
        eval = double(f(x(j,i),a(i)));
        format long
        error = abs(eval);
        
        message1 = strcat('x = ',num2str(x(j,i)),' f(x) = ',num2str(eval));
        disp(message1)
    end
    root(i) = x(j,i);
    errorTracker(i) = error;
    numIterations(i) = j;
    message = strcat(num2str(root(i))," is a root!");
    disp(message)
    message4 = ('--------------------------------');
    disp(message4)
end
end % End of Bisection Method

