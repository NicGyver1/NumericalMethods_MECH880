% Author: Nicholas Piercy
% Date : 11/17/2021
% Numerical Methods, Homework #6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%% Problem 1
interval = [-1 2]; % Interval to approximate over
n = [0 1 2 3 4];
for i = 1:length(n)
    myfunc = @(x) x^n(i); % Function to approximate
    I = SingleSimpsons(interval,myfunc);
    disp(I)
end

%% Problem 2
interval = [1 4];
myfunc = @(x) exp(-x)*sin(x);

Exact = (-exp(-4)/2)*(sin(4)+cos(4)) - (-exp(-1)/2)*(sin(1)+cos(1));

numintervals = [1 2 4 8 16 32 64 128];
disp('Simpsons Method')
for i = 1:length(numintervals)
    I_simp = Composite_Simpsons(interval,myfunc,numintervals(i));
    rel_error_simp(i) = (abs(Exact-I_simp)/Exact);
    message = strcat('h = ',num2str(numintervals(i)),', ',' function value = '...
        ,num2str(I_simp),' Rel Error = ',num2str(rel_error_simp(i)));
    disp(message)
end
disp('Midpoint Method')
for i = 1:length(numintervals)
    I_midpt = Composite_Midpoint(interval,myfunc,numintervals(i));
    rel_error_midpt(i) = (abs(Exact-I_midpt)/Exact);
    message = strcat('h = ',num2str(numintervals(i)),', ',' function value = ',...
        num2str(I_midpt),' Rel Error = ',num2str(rel_error_midpt(i)));
    disp(message)
end
disp('Trapezoidal Method')
for i = 1:length(numintervals)
    I_trap = Composite_Trapezoidal(interval,myfunc,numintervals(i));
    rel_error_trap(i) = (abs(Exact-I_trap)/Exact);
    message = strcat('h = ',num2str(numintervals(i)),', ',' function value = ',...
        num2str(I_trap),' Rel Error = ',num2str(rel_error_trap(i)));
    disp(message)
end

plot(log(numintervals),log(rel_error_simp),log(numintervals),log(rel_error_midpt),...
    log(numintervals),log(rel_error_trap))
legend('Simspons','Midpoint','Trapezoidal')
title('Loglog plot for error convergence','interpreter','latex');
ylabel('log(Error)','interpreter','latex');
xlabel('log(h)','interpreter','latex')
%% Problem 3

%% Problem 4
clc

exact = pi;
interval = [0 1];
myfunc4 = @(x) 4./(1+x.^2);
matlab_exact = integral(myfunc4,interval(1),interval(2));

xi1 = -sqrt(3/5); xi2 = 0; xi3 = sqrt(3/5); % gaussian location points
w1 = 5/9; w2 = 8/9; w3 = 5/9; % corresponding gaussian weights

a = interval(1);
b = interval(2);
numintervals = [1 2 4 8 16];
message1 = strcat('Num Intervals,',' Approximation,',' Matlab Numerical Integration',' Exact Error');
disp(message1)
for i = 1:length(numintervals)
    h(i) = (interval(2)-interval(1))/numintervals(i);
    x = linspace(interval(1),interval(2),(numintervals(i)+1));
    I = 0;
    for j = 1:numintervals(i)
        a = x(j);
        b = x(j+1);
        n = (a+b)/2;
        m = (b-a)/2;
        I = I + m*w1*my_fp4(m*xi1 + n) + m*w2*my_fp4(m*xi2 + n) + m*w3*my_fp4(m*xi3 + n);
    end
    Integral(i) = I;
    Error(i) = abs(exact - Integral(i));
    message = [numintervals(i), Integral(i), matlab_exact, Error(i)];
    disp(message);
end
plot(log(h),log(Error))
legend('3-noded Gaussian Quad')
title('Gaussian Quadrature Error','interpreter','latex');
ylabel('log(Error)','interpreter','latex');
xlabel('log(h)','interpreter','latex')


%% Problem 5
x=0.66; 
% the exact value of the derivative at x
exact=(exp(x)*(x-3))./(x-2)^2;
h=0.01; % the step size
for i=1:3
fw=(my_f(x+h)-my_f(x))/h; % forward difference approximation
bw=(my_f(x)-my_f(x-h))/h; % backward difference approximation
cd=(my_f(x+h)-my_f(x-h))/(2*h); % central difference approximation
errorfw=exact-fw;
errorbw=exact-bw;
errorcd=exact-cd;
s=[h, fw, bw, cd, errorfw, errorbw, errorcd]; % print the step size, the approx. value, and the absolute error
disp(s);
h=h/10;
end

%% Problem 6
x = [1 pi];
flag = 1;
for i = 1:length(x)
    exact = (-4*sin(x(i)))-(2*x(i)*cos(x(i)));
    h = 0.1;
    for j = 1:6
        cd_dp = (my_fp6(x(i)+h) - 2*my_fp6(x(i)) + my_fp6(x(i)-h))/(h^2);
        error_cd_dp(flag) = abs(exact-cd_dp);
        s = [h, cd_dp, error_cd_dp(flag)];
        disp(s);
        
        htracker(flag) = h;
        flag = flag+1;
        h = h/10;
    end
end
 plot(log(htracker(1:6)),log(error_cd_dp(1:6)))
 hold on   
 plot(log(htracker(7:12)),log(error_cd_dp(7:12)))
 legend('at x = 1','at x = $\pi$','interpreter','latex');
 title('Central Difference Error for Approximating $f(x)=2xcos(x)$','interpreter','latex');
 ylabel('log(Error)','interpreter','latex');
 xlabel('log(h)','interpreter','latex')

%% Functions
function [y] = my_f(x);
y=exp(x)./(x-2);
end
function [yyy] = my_fp4(x);
yyy = 4/(1+x^2);
end
function [yy] = my_fp6(x);
yy = 2*x*cos(x);
end

