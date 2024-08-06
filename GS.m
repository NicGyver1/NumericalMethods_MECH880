% Author: Nicholas Piercy
% Date : 9/30/2021
% Gauss-Seidel Iterative Method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,k,epsilon] = GS(A,b,x_guess,hard_stop,epsilon_stop,normNum)
format long
n = length(A(1,:));
if x_guess == 0
    x_guess = zeros(n,1);
end
x = x_guess;
for k = 1:hard_stop
    for i = 1:n
        x(i) = (b(i) - (dot(A(i,1:i-1),x(1:i-1))...
            + dot(A(i,i+1:n),x(i+1:n))))/A(i,i);
    end    
    epsilon = norm((A*x-b),normNum); % norm of the relative error
    if epsilon < epsilon_stop
        break
    end
end
message1 = char({['Gauss Seidel Iterative Method'];
            ['Stopping Criterion Tolerance = ',num2str(epsilon_stop)];
            ['Number of Iterations To Converge  = ', num2str(k)];
            ['Norm of Residual = ', num2str(epsilon)]});
if k == hard_stop
    message1 = char({[message1];['Did not converge, try another guess!']});
    disp(message1)
else 
    disp(message1) % print to screen the number of iterations to converge
end
for i = 1:length(x)
    message2 = char({['x(',num2str(i),') =  ', num2str(x(i))]});
    disp(message2)
end
