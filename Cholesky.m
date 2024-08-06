% Author: Nicholas Piercy
% Date : 9/18/2021
% Cholesky Factorization Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H,HT] = Cholesky(A)
n = length(A(:,1));
for k = 1:n
    A(k,k) = sqrt(A(k,k)); %Since diagonal terms will be multiplied together
    for i = k+1:n
        A(k,i) = A(i,k)/A(k,k); % Since A(k,k) term will be multiplied back in
    end % First row is now factored 
    for j = k+1:n
        for i = j:n % finding Lower diagonal
            A(j,i) = A(i,j)-A(k,i)*A(k,j);
        end
    end
    A(k+1:n,k) = 0;
end
HT = A;
H = HT'; 
end % end of Cholesky function
