% Author: Nicholas Piercy
% Date : 9/18/2021
% Gauss Jordan method utilizing Identity matrix to find inverse.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
% A     --> NonSingular, Invertible Matrix
% b     --> RHS Vector(or matrix)

%OUTPUTS
% D     --> Diagonal Matrix 
% b_new --> Corresponding RHS Vector(or matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Testing
% clear all
% A = [2 2 4; 1 2 6; 1 2 8];
% b = [2;5;6];

% [Inverse,b_new,check] = GaussJordan_Inverse(A,b);
%%

function [A_Inverse,Check,varargout] = GaussJordanInverseV3(varargin)
if nargin == 1
%     vargout = 2;
    A = varargin{1};
    [A_Inverse,Check] = GaussJordanInverse_withoutRHSV3(A);
%     vargout(1) = Inverse;
    varargout{1} = [];
elseif nargin == 2
%     vargout = 3;
    A = varargin{1};
    b = varargin{2};
    [A_Inverse,b_new,Check] = GaussJordanInverse_withRHSV3(A,b);
%     vargout(1) = Inverse;
    varargout{1} = b_new;
%     vargout = Check;
end
end % end of "GaussJordan" function

%% If only Matrix A is given
function [Inverse,Check] = GaussJordanInverse_withoutRHSV3(A)

n = length(A(:,1)); % number of rows in matrix
Inverse = eye(n);
Check = inv(A);
AI = [A,Inverse];


if nargin == 1
for i = 1:(n-1)
    [~,index] = max(abs(A(i:n,i))); % find row index for largest value (only considering rows below previous)
    index = index + i - 1; % focuses on the submatrix which hasnt seen pivots yet
    %perform partial pivoting
    if index ~= i
        %swapping pivot rows
        AI_move = AI(i,i:end); % temporarily storing old pivot row
        AI(i,i:end) = AI(index,i:end); % replacing old pivot row with new
        AI(index,i:end) = AI_move; % filling new pivot row's old spot with old pivot row
    end
    %compute multipliers
    for j = i+1:n
        m_L(j) = AI(j,i)/AI(i,i);
    end
    for j = i+1:n % loop over rows
        for k = 1:length(AI(1,:)) % loop over columns
            AI(j,k) = AI(j,k) - m_L(j)*AI(i,k);
        end
    end
    %Compute multipliers for upper diagonal elimination
    for j = 1:i
        m_U(j) = AI(j,i+1)/AI(i+1,i+1);
    end
    for j = 1:i % loop over rows
        for k = i+1:length(AI(1,:)) % loop over columns
            AI(j,k) = AI(j,k) - m_U(j)*AI(i+1,k);
%             Inverse(j,k) = Inverse(j,k) - m_U(j)*Inverse(i+1,k);
        end
    end 
end
for i = 1:n
    AI(i,:) = AI(i,:)/AI(i,i); %Divide  [A|I] by pivot magnitude
end
Inverse = AI(:,n+1:end);
end
end



%% If Matrix A AND Vector/Matrix b is given
function [Inverse,b_new,Check] = GaussJordanInverse_withRHSV3(A,b)
n = length(A(:,1)); % number of rows in matrix
Inverse = eye(n);
Check = inv(A);
AI = [A,Inverse];
b_size = size(b); % check if b is a vector or matrix


for i = 1:(n-1)
    [~,index] = max(abs(A(i:n,i))); % find row index for largest value (only considering rows below previous)
    index = index + i - 1; % focuses on the submatrix which hasnt seen pivots yet
    
    %perform partial pivoting
    if index ~= i
        %swapping pivot rows
        AI_move = AI(i,i:end); % temporarily storing old pivot row
        AI(i,i:end) = AI(index,i:end); % replacing old pivot row with new
        AI(index,i:end) = AI_move; % filling new pivot row's old spot with old pivot row
        
        %swapping RHS vector(or matrix) corresponding to swap in LHS
        b_move = b(i,i:end);
        b(i,i:end) = b(index,i:end);
        b(index,i:end) = b_move;
    end
    
    %compute multipliers
    for j = i+1:n
        m_L(j) = AI(j,i)/AI(i,i);
    end
    
    for j = i+1:n % loop over rows
        for k = 1:length(AI(1,:)) % loop over columns
            AI(j,k) = AI(j,k) - m_L(j)*AI(i,k);
%             Inverse(j,k) = Inverse(j,k) - m_L(j)*Inverse(i,k);
            if b_size(2) ~=1 % account for RHS as matrix
                b(j,k) = b(j,k) - m_L(j)*b(i,k);
            end
        end
        if b_size(2) ==1
            b(j) = b(j) - m_L(j)*b(i);
        end
    end
    
    %Compute multipliers for upper diagonal elimination
    for j = 1:i
        m_U(j) = AI(j,i+1)/AI(i+1,i+1);
    end

    for j = 1:i % loop over rows
        for k = i+1:length(AI(1,:)) % loop over columns
            AI(j,k) = AI(j,k) - m_U(j)*AI(i+1,k);
%             Inverse(j,k) = Inverse(j,k) - m_U(j)*Inverse(i+1,k);
            if b_size(2) ~=1 % account for RHS as matrix
                b(j,k) = b(j,k) - m_U(j)*b(i+1,k);
            end
        end
        if b_size(2) ==1
            b(j) = b(j) - m_U(j)*b(i+1);
        end
    end 
end
for i = 1:n
    AI(i,:) = AI(i,:)/AI(i,i); %Divide  [A|I] by pivot magnitude
    b(i,:) = b(i,:)/AI(i,i); %Divide RHS by pivot magnitude
end
    
Inverse = AI(:,n+1:end);
b_new = b;
end


   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    