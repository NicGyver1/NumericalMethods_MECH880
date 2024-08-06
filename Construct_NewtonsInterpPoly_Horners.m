% Author: Nicholas Piercy
% Date : 11/1/2021
% Contrsuct Newtons Polynomial using Horners Nested Evaluation Scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
% gamma    --> Coefficients from Divided Difference Table
% x_values --> x values of data input 

%OUTPUTS
% Pn    --> Interpolating Polynomial 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Horners Nestd Evaluation Scheme
function PnHorners = Construct_NewtonsInterpPoly_Horners(gamma,x_values)
% Pn(x) = x(1) + (x - x(1))*(x(2) + (x - x(2)))*(x(3) + (x - x(3)))...
% So

for i = 1:n-1
    if i == 1
        prod = (x - x_vals(i));
    elseif i < n-1
        prod = prod*(gamma(i) + (x - x_vals(i)));
    else
        prod = prod*(gamma(i) + gamma(i+1)*(x - x_vals(i)));
    end
end
PnHorners(x) = gamma(1) + prod