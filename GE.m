% Author: Nicholas Piercy
% Date : 9/14/2021
% Gaussian Elimination Without Partial Pivoting (GE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
% A     --> NonSingular, Invertible Matrix
% b     --> RHS Vector(or matrix)

%OUTPUTS
% U     --> Upper Traingular Matrix 
% b_new --> Corresponding RHS Vector(or matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,b_new] = GE(A,b)
n = length(A); % number of rows in matrix
b_size = size(b); % check if b is a vector or matrix
for i = 1:(n-1)
%     [~,index] = max(abs(A(i:end,i))); % find row index for largest value (only considering rows below previous)
%     index = index + i - 1; % focuses on the submatrix which hasnt seen pivots yet
    
    %compute multipliers
    for j = i+1:n
        m(j) = A(j,i)/A(i,i);
    end
    for j = i+1:n % loop over rows
        for k = 1:n % loop over columns
            A(j,k) = A(j,k) - m(j)*A(i,k);
            if b_size(2) ~=1 % account for RHS as matrix
                b(j,k) = b(j,k) - m(j)*b(i,k);
            end
        end
        if b_size(2) ==1
            b(j) = b(j) - m(j)*b(i);
        end
    end
end
U = A;
b_new = b;
end % end of "GE" function

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    