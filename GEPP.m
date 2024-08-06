% Author: Nicholas Piercy
% Date : 9/14/2021
% Gaussian Elimination With Partial Pivoting (GEPP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
% A     --> NonSingular, Invertible Matrix
% b     --> RHS Vector(or matrix)

%OUTPUTS
% U     --> Upper Traingular Matrix 
% b_new --> Corresponding RHS Vector(or matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,b_new] = GEPP(A,b)

n = length(A); % number of rows in matrix
b_size = size(b); % check if b is a vector or matrix
for i = 1:(n-1)
    [~,index] = max((abs(A(i:n,i)))); % find row index for largest value (only considering rows below previous)
    currentRow = i-1;
    index = index + currentRow; % focuses on the submatrix which hasnt seen pivots yet
    
    %perform partial pivoting
    if index ~= i
        %swapping pivot rows
        A_move = A(i,i:n); % temporarily storing old pivot row
        A(i,i:n) = A(index,i:n); % replacing old pivot row with new
        A(index,i:n) = A_move; % filling new pivot row's old spot with old pivot row
        
        
        
        
        %swapping RHS vector(or matrix) corresponding to swap in LHS
        b_move = b(i,:);
        b(i,:) = b(index,:);
        b(index,:) = b_move;
    end
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
end % end of "GEPP" function

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    