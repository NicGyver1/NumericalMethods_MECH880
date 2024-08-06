% Author: Nicholas Piercy
% Date : 11/1/2021
% Contrsuct Newtons Polynomial.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
% gamma    --> Coefficients from Divided Difference Table
% x_values --> x values of data input 

%OUTPUTS
% Pn    --> Interpolating Polynomial 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Evaluate Newtons Interpolating polynomial
function Pn = Construct_NewtonsInterpPoly(gamma,x_values)
syms x Pn(x) unknowns
n = length(x_values);
unknowns(1) = 1;
for i = 2:n
    unknowns(i) = unknowns(i-1)*(x - x_values(i-1));
end
Pn(x) = dot(gamma,unknowns);
end % end of Construct_NewtonsInterpPoly function