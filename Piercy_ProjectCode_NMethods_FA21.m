% Author: Nicholas Piercy
% Date : 12/12/2021
% Numerical Methods, Project "Computing Eigenvalues and Eigenvectors"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given Matrix 
A = [5 4 2 1;
    0 1 -1 -1;
    -1 -1 3 0;
    1 1 -1 2];
% Initial Eigenvector Guess
x0 = [1;0;1;0]; 
I = eye(length(A(:,1)));
% Rayleigh Quotient for Initial Eigenvalue Guess
q = (x0'*(A*x0))/(x0'*x0); 
for i = 1:50
    AqI = A - q*I; %Setup Inverse Power Method
    y = AqI\x0; % solve linear system
    [~,index] = max(abs(y));
    y_pm = y(index); % Find Index of Largest Magnitude Entry
    x_m = y/y_pm; %Normalized Eigenvector Approximation
    lambda = q + 1/y_pm; %Eigenvalue Approximation
    error = abs(x0 - x_m);
    %stopping criterion
    if error < 1e-5; break; end
    %Set up next iteration
    x0 = x_m;
    %New Eigenvalue Guess
    q = 1/(lambda-q); % For Smallest Eigenvalue
%     q = (lambda-q); % For Largest Eigenvalue
end
IterationsToConverge = i;
message1 = ("The Eigenvalue closest to q is:");
message2 = ("The corresponding eigenvector is:");
message3 = ("Number of iterations to converge:");
disp(message1)
disp(num2str(lambda))
disp(message2)
disp(num2str(x_m(1)))
disp(num2str(x_m(2)))
disp(num2str(x_m(3)))
disp(num2str(x_m(4)))
disp(message3)
disp(IterationsToConverge)
[v,d] = eig(A); % checking solution
disp(v)
disp(d)

