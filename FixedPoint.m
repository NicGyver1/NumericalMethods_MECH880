% Author: Nicholas Piercy
% Date : 10/19/2021
% Fixed-Point method for finding a root of a continuous function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% interval  -> Interval which contains a root, [a,b]
% f         -> Function
% tol       -> Error tolerance for convergence
function [x,error,j] = FixedPoint(guess,g,f,tol)
error = 1; % Initial value of error (dummy)
eval_previous = 0; % for Problem 8_HW4 numerical methods UNL 10/19/2021
j = 1;
while error>tol
    
    if j ==1
        x(j) = g(guess);
    else
        x(j) = g(x(j-1));
    end
    eval = double(f(x(j)));
    
    error_P8_HW4 = abs(eval-eval_previous);
    error = error_P8_HW4;
    
    if error_P8_HW4<tol
        break
    end
    eval_previous = eval;
    
    %Print out iterate information
    message1 = strcat('x = ',num2str(x(j)),' f(x) = ',num2str(eval)...
        ,'    |x(k+1)-x(k)| = ',num2str(error));
    disp(message1)
    j = j + 1;
end
message = strcat(num2str(x(j))," is a root!");
disp(message)
message2 = strcat('Iterations to converge = ',num2str(j));
message3 = strcat('Error at termination = ',num2str(error_P8_HW4));
disp(message2)
disp(message3)
end % End of Bisection Method

