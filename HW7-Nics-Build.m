% Author: Nicholas Piercy
% Date : 11/17/2021
% Numerical Methods, Homework #6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% Problem 1
clear all
close all
clc

y_previous = 2; % y(0) = 2 but with index starting at 1...
interval = [0 1];
h = [10^(-1) 10^(-2) 10^(-3) 10^(-4)];
for i = 1:length(h)
    timesteps = (interval(2)-interval(1))/h(i);
    
    y_next = zeros(timesteps,1);
    for t = 0:1:timesteps-1
        func1 = t + y_previous - t*y_previous^3;
        y_next(t+1) = ForwardEuler(func1,h(i),y_previous);
        y_previous = y_next(t+1);
    end
    t = linspace(interval(1),interval(2),timesteps); % points for plotting
    plot(t,y_next)
    xlim = [interval(1) interval(2)];
    hold on
end
title("Forward Euler for $y'(t) = t + y(t) + ty^3(t)$",'fontsize',15,'interpreter','latex')
xlabel('$t$','fontsize',15,'interpreter','latex')
ylabel("$y'(t)$",'fontsize',15,'interpreter','latex')
legend('$h=10^{-1}$','$h=10^{-2}$','$h=10^{-3}$','$h=10^{-4}$','fontsize',15,'interpreter','latex')
% unstable, especially for small step sizes near $t \approx 1$;  
%  --> This is expected since this is an explicit method... => less calcs but less precise
    

%% Problem 2
clear all
close all
clc

y_previous = 0; % 0 initial velocity at the jump out of the airplane!
h = 0.2; % time step size
interval = [0:h:3];% time to calculate
y_next = zeros((length(interval)),1);
for i = 2:length(interval)
    func2 = myfunc2(y_previous);
    Predictor = ForwardEuler(func2,h,y_previous);
    Corrector = y_previous + (h/2)*(func2 + myfunc2(Predictor));
%     func2 = g*(1-(y_previous/80)^(3/2));
    y_next(i) = Corrector;
    y_previous = y_next(i);
    
end
VelocityAt3Seconds = y_next(end);
plot(interval,y_next)
hold on

y0 = 0;
tspan = [0 3];
g = 9.8;
exactSoln = ode45(@(t,y) (g*(1-((y/80)^(3/2)))), tspan, y0);
plot(exactSoln.x,exactSoln.y,'o')

%% Problem 4
clear all
close all
clc

y1_0 = 100;
y2_0 = 10;
initialConditon = [y1_0;y2_0];

% f1 = y1*(1-0.1*y2);
% f2 = y2*(-0.5 + 0.2*y1);

% func = [f1;f2];
tspan = [0 25];

[t,y] = ode45(@predprey4,tspan,initialConditon);
figure
plot(t,y(:,1),t,y(:,2))
legend('$t$ vs. $y_1$,$y_2$','interpreter','latex')
figure
plot(y(:,1),y(:,2))
legend('$y_1$ vs. $y_2$','interpreter','latex')

function dydt = predprey4(t,y)
alpha1 = 1;
alpha2 = 0.5;
beta1 = 0.1;
beta2 = 0.02;

dydt = zeros(2,1);


dydt(1) = y(1)*(alpha1-beta1*y(2));
dydt(2) = y(2)*(-alpha2 + beta2*y(1));
% dydt = [y(1)*(alpha1-beta1*y(2));y(2)*(-alpha2 + beta2*y(1))];
end


%% Functions
function [y_next] = ForwardEuler(func,h,y_previous)
y_next = y_previous + h*func;
end

function [y_next] = myfunc2(y_previous)
g = 9.8; %m/(s^2) acceleration due to gravity
y_next = g*(1-(y_previous/80)^(3/2));
end