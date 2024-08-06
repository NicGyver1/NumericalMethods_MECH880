% Author: Nicholas Piercy
% Date : 9/30/2021
% Jacobi Iterative Method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_new,k,epsilon] = Jacobi(A,b,x_guess,hard_stop,epsilon_stop,normNum)
format long
n = length(A(1,:));
if x_guess == 0
   x_guess = zeros(n,1);
end
x_new = zeros(n,1);
for k = 1:hard_stop
    for i = 1:n
        x_new(i) = (b(i) - (dot(A(i,1:i-1),x_guess(1:i-1))...
            + dot(A(i,i+1:n),x_guess(i+1:n))))/A(i,i);
    end
    epsilon = norm((A*x_new-b),normNum);
    if epsilon < epsilon_stop
        break
    end
    x_guess = x_new;
end
message1 = char({['Jacobi Iterative Method'];
            ['Stopping Criterion Tolerance = ',num2str(epsilon_stop)];
            ['Number of Iterations To Converge  = ', num2str(k)];
            ['Norm of Residual = ', num2str(epsilon)]});
if k == hard_stop
    message1 = char({[message1];['Did not converge, try another guess!']});
    disp(message1)
else 
    disp(message1) % print to screen the number of iterations to converge
end
for i = 1:length(x_new)
    message2 = char({['x(',num2str(i),') =  ', num2str(x_new(i))]});
    disp(message2)
end
