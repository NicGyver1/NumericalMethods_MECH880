% Author: Nicholas Piercy
% Date : 10/18/2021
% Numerical Methods, Homework #4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 2
clear all
clc
syms t
X = [0,1,3];
Y = [1,2,3];

A = [1 0;
    1 1;
    1 3];
b = [1; 2; 3];
[M,d] = NormalEqs(A,b);
[H,HT] = Cholesky(M);
y = ForwardSub(H,d);
x = BackwardSub(HT,y);
x_checker = M\d;
disp(x)

func = x(1) + x(2)*t;
domain = [-1 4];
fplot(func,domain)
hold on
plot(X,Y,'o') 
axis equal
grid on
xlim([-1 4]); ylim([0 4])
xlabel('$x$','interpreter','latex'); 
ylabel('$f(x)$','interpreter','latex');
legend('Linear Fit','fontsize',15,'interpreter','latex','linewidth',2)
    
%% Problem 4
clear all; clc;
syms x xp

t = [0, 1, 2, 3, 4, 5];
y = [1; 2.7; 5.8; 6.6; 7.5; 9.9];
A = [1; 1; 1; 1; 1; 1];
for n = 1:5
    for i = 1: length(t)
        A(i,n+1) = (t(i))^n;
    end
    [M,d] = NormalEqs(A,y);
    [U,d_new] = GE(M,d);
    x_gepp{n} = BackwardSub(U,d_new);
   
    xp(1) = 1;
    for k = 1:n
        xp(k+1) = x^k;
    end
    
    func = dot(x_gepp{n},xp);
    xs = [min(t)-1, max(t)+1];
    if n==1; width = 4; % Show Linear fit is best  
    else; width = 1;
    end
    fplot(func,xs,'Linewidth',width)
    hold on
end
plot(t,y,'ro','linewidth',5,'markersize',5)
grid on
xlim([-1 6]); ylim([0 11])
xlabel('$x$','interpreter','latex'); 
ylabel('$f(x)$','interpreter','latex');
legend('Linear Fit','Quadratic Fit','Cubic Fit'...
    ,'Quartic Fit','Quintic Fit','Data Points',...
    'fontsize',15,'interpreter','latex','linewidth',2)

 

%% Problem 6
clear all
clc
syms x 

f_6 = @(x)(sqrt(x)-1.1); % function
interval_6 = [0,2];  % interval [a,b]
tol_6 = 10^(-8);
[root_6,error_6,numIterations_6] = bisection(interval_6,f_6,tol_6);

%% Problem 7
clear all
clc
syms f(x,a)
f(x,a) = (x^3 - a);

a_7 = [0, 2, 10];
guess_7 = [.025, 1.25, 2.2];

tol_7 = 10^(-8);
[root_7,error_7,numIterations_7,x] = NewtonsMethod(guess_7,a_7,f,tol_7);

%% Problem 8
clear all
clc
syms x func(x) g

f = @(x) (x + log(x));
g = @(x) (exp(-x));
domain = [0.1 1];
fplot(f,domain,'linewidth',2);
grid on
xlabel('$x$','interpreter','latex'); 
ylabel('$f(x)$','interpreter','latex');
legend('$f(x)=x+ln(x)$','interpreter'...
    ,'latex','fontsize',15,'linewidth',2)

interval = [0.5 0.6];
tol = 10^(-10);

[root_bisection,error_bisection,iterations_bisection] = bisection(interval,f,tol);

guess_FixedPoint = 0.5;
[root_FixedPoint,error_FixedPoint,iterations_FixedPoint] = FixedPoint(guess_FixedPoint,g,f,tol);

func(x) = (x + log(x));


newton_guess = (0.5);
[root_newton,error_newton,iterations_newton,x_newton] = NewtonsMethod_oneVar(newton_guess,func,tol);

secant_guess = [0.5 0.6];
[root_secant,error_Secant,iterations_secant,x_secant] = SecantMethod(secant_guess,func,tol);



