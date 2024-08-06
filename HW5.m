% Author: Nicholas Piercy
% Date : 11/1/2021
% Numerical Methods Homework #5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

%% Problem 2
% Compute Newton Divided Difference Table and get Coefficients
[gamma,DividedDifference_Table,n] = Newton_DivDiff(x_vals_2,y_vals_2);

% Construct Newton Basis Polynomial
Pn = Construct_NewtonsInterpPoly(gamma,x_vals_2); 

% Evaluate the Polynomial at the point "x_evaluate"
% Pn(x_evaluate); 

 % Compute new Coefficients for added point
[gamma_new,DivDiff_new,x_vals_added] = Newton_DivDiff_PointAdded(x_new...
    ,y_new,x_vals_2,DividedDifference_Table,gamma);

% x_new = 7;
% y_new = 48;
% [gamma_new,DivDiff_new,x_vals_added] = Newton_DivDiff_PointAdded(x_new,y_new,x_vals,DividedDifference_Table,gamma);
% Pn_added = Construct_NewtonsInterpPoly(gamma_new,x_vals_added);
% figure
% fplot(@(x) Pn_added(x),[0 8])
%% Problem 3
x_vals_3 = [0 0.5 1 6 7.0 9.0];
y_vals_3 = [0 1.6 2.0 2.0 1.5 0.0];
[gamma,DividedDifference_Table,n] = Newton_DivDiff(x_vals_3,y_vals_3);

Pn = Construct_NewtonsInterpPoly(gamma,x_vals_3)
figure
fplot(@(x) Pn(x),[0 9])
hold on

xqp3 = (0:0.1:9);
s_cubicspline_3 = interpn(x_vals_3,y_vals_3,xqp3,'pchip');
s_piecewiselinear_3 = interpn(x_vals_3,y_vals_3,xqp3,'linear');
plot(x_vals_3,y_vals_3,'o',xqp3,s_piecewiselinear_3,'-.',xqp3,s_cubicspline_3,'-')
legend('Newtons Interpolation','Data Points','Linear Interpolant','Cubic Spline','location','best')

%% Problem 4
x_vals_4 = [0 10 20 30 40 50 60 70 80];
x_vals_4 = 1900+x_vals_4;
y_vals_4 = [76212168  92228496 106021537 123202624 132164569 151325798 179323175 203302031 226542199];

xqp4 = (1900:1:1980); % Domain for plotting

% Splines
s_pchipspline_4 = interpn(x_vals_4,y_vals_4,xqp4,'pchip');
s_cubicspline_4 = interpn(x_vals_4,y_vals_4,xqp4,'cubic');


% Newtons Interp Poly
[gamma4,DividedDifference_Table4,n] = Newton_DivDiff(x_vals_4,y_vals_4);
Pn_4 = Construct_NewtonsInterpPoly(gamma4,x_vals_4);

fplot(@(x) Pn_4(x),[1900 1990])
hold on
% Newtons interp with added point from spline
x_new4 = 1990;
y_new4 = spline(x_vals_4,y_vals_4,1990);
[gamma_new4,DivDiff_new,x_vals_added4] = Newton_DivDiff_PointAdded(x_new4,y_new4,x_vals_4,DividedDifference_Table4,gamma4);
Pn_added4 = Construct_NewtonsInterpPoly(gamma_new4,x_vals_added4);
fplot(@(x) Pn_added4(x),[1900 1990])
hold on

% Plot splines
plot(x_vals_4,y_vals_4,'o','linewidth',2)
plot(xqp4,s_pchipspline_4,'-.','linewidth',2)
plot(xqp4,s_cubicspline_4,'-','linewidth',2)
hold on

legend('Newtons Pn8','Newtons Pn9','data points','Hermite Spline','cubic spline')





