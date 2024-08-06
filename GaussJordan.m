% Author: Nicholas Piercy
% Date : 9/14/2021
% Gauss Jordan method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
% A     --> NonSingular, Invertible Matrix
% b     --> RHS Vector(or matrix)

%OUTPUTS
% D     --> Diagonal Matrix 
% b_new --> Corresponding RHS Vector(or matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D,b_new] = GaussJordan(A,b)
n = length(A); % number of rows in matrix
b_size = size(b); % check if b is a vector or matrix
for i = 1:(n-1)
    [~,index] = max(abs(A(i:end,i))); % find row index for largest value (only considering rows below previous)
    index = index + i - 1; % focuses on the submatrix which hasnt seen pivots yet
    %perform partial pivoting
    if index ~= i
        %swapping pivot rows
        A_move = A(i,i:end); % temporarily storing old pivot row
        A(i,i:end) = A(index,i:end); % replacing old pivot row with new
        A(index,i:end) = A_move; % filling new pivot row's old spot with old pivot row
        %swapping RHS vector(or matrix) corresponding to swap in LHS
        b_move = b(i,i:end);
        b(i,i:end) = b(index,i:end);
        b(index,i:end) = b_move;
    end
    %compute multipliers
    for j = i+1:n
        m_L(j) = A(j,i)/A(i,i);
    end
    for j = i+1:n % loop over rows
        for k = 1:n % loop over columns
            A(j,k) = A(j,k) - m_L(j)*A(i,k);
            if b_size(2) ~=1 % account for RHS as matrix
                b(j,k) = b(j,k) - m_L(j)*b(i,k);
            end; end;
        if b_size(2) ==1
            b(j) = b(j) - m_L(j)*b(i);
        end; end;
    %Compute multipliers for upper diagonal elimination
    for j = 1:i
        m_U(j) = A(j,i+1)/A(i+1,i+1);
    end
    for j = 1:i % loop over rows
        for k = i+1:n % loop over columns
            A(j,k) = A(j,k) - m_U(j)*A(i+1,k);
            if b_size(2) ~=1 % account for RHS as matrix
                b(j,k) = b(j,k) - m_U(j)*b(i+1,k);
            end; end;
        if b_size(2) ==1
            b(j) = b(j) - m_U(j)*b(i+1);
        end; end; end;
D = A;
b_new = b;
end % end of "GaussJordan" function

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    