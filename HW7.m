% Author: Nicholas Piercy
% Date : 12/2/2021
% Numerical Methods, Homework #7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% Problem 1
y0 = 2; % Initial Condition
tspan = [0 1];
h = [10^(-1) 10^(-2) 10^(-3) 10^(-4)];
for i = 1:length(h)
    timesteps = (tspan(2)-tspan(1))/h(i);
    t{i} = linspace(tspan(1),tspan(2),timesteps); % points for plotting
    y{i} = zeros(timesteps,1);
    y{i}(1) = y0;
    for j = 2:timesteps
        func1 = t{i}(j) + y{i}(j-1) - t{i}(j)*y{i}(j-1)^3;
        y{i}(j) = ForwardEuler(func1,h(i),y{i}(j-1));
    end
end
plot(t{1},y{1},t{2},y{2},t{3},y{3},'-',t{4},y{4},'linewidth',2)
hold on
[t_matlab,y_matlab] = ode45(@(t,y) t+y-t*y^3,tspan,y0);
plot(t_matlab,y_matlab,'--d','linewidth',2)
xlim = [tspan(1) tspan(2)];
title("Forward Euler for $y'(t) = t + y(t) + ty^3(t)$",'fontsize',15,...
    'interpreter','latex')
xlabel('$t$','fontsize',15,'interpreter','latex')
ylabel("$y'(t)$",'fontsize',15,'interpreter','latex')
legend('$h=10^{-1}$','$h=10^{-2}$','$h=10^{-3}$','$h=10^{-4}$',...
    '$y_{ode45}$','fontsize',15,'interpreter','latex')
%% Problem 2
y_previous = 0; % 0 initial velocity at the jump out of the airplane!
h = 0.2; % time step size
tspan = [0:h:3];% time to calculate
y = zeros((length(tspan)),1);
for i = 2:length(tspan)
    func2 = myfunc2(y_previous);
    Predictor = ForwardEuler(func2,h,y_previous);
    Corrector = y_previous + (h/2)*(func2 + myfunc2(Predictor));
    y(i) = Corrector;
    y_previous = y(i);
end
VelocityAt3Seconds = y(end);
plot(tspan,y)
hold on
y0 = 0;
tspan = [0 3];
g = 21.937; % [mph/s]
exactSoln = ode45(@(t,y) (g*(1-((y/80)^(3/2)))), tspan, y0);
plot(exactSoln.x,exactSoln.y,'o')
legend('Predictor-Corrector','$y_{ode45}$','fontsize',15,'interpreter','latex')
xlabel('$t$','fontsize',15,'interpreter','latex')
ylabel("$Velocity$",'fontsize',15,'interpreter','latex')
%% Problem 4
y1_0 = 25;
y2_0 = 10;
initialConditon = [y1_0;y2_0];
tspan = [0 25];
[t,y] = ode45(@predprey4,tspan,initialConditon);
figure
plot(t,y(:,1),t,y(:,2)) % Solution
xlabel('t','interpreter','latex')
ylabel('$y_i(t)$','interpreter','latex')
legend('Prey','Predator')
title('Function Plot, $y_1(0)=25,y_2(0)=10$','interpreter','latex')
figure
plot(y(:,1),y(:,2)) % Phase Portrait
xlabel('$y_1(t)$','interpreter','latex')
ylabel('$y_2(t)$','interpreter','latex')
legend('$y_1$ vs. $y_2$','interpreter','latex')
title('Phase Portrait, $y_1(0)=25,y_2(0)=10$','interpreter','latex')
%% Problem 4 - part 2
alpha1 = 1; alpha2 = 0.5; beta1 = 0.1; beta2 = 0.02;
h = [10^(-2), 10^(-4)]; % time step size
tspan = [0 25];
for j = 1:length(h)
timesteps = (tspan(2)-tspan(1))/h(j);
t{j} = linspace(tspan(1),tspan(2),timesteps); % points for plotting
y1{j} = zeros(timesteps,1); y2{j} = zeros(timesteps,1);
y1{j}(1) = 100; y2{j}(1) = 10; % initial conditions
    for i = 2:timesteps
        func4_1 = y1{j}(i-1)*(alpha1-beta1*y2{j}(i-1));
        func4_2 = y2{j}(i-1)*(-alpha2 + beta2*y1{j}(i-1));
        y1{j}(i) = ForwardEuler(func4_1,h(j),y1{j}(i-1));
        y2{j}(i) = ForwardEuler(func4_2,h(j),y2{j}(i-1));
    end
end
figure
plot(t{1},y1{1},t{1},y2{1},t{2},y1{2},t{2},y2{2})
xlabel('t','interpreter','latex')
ylabel('$y_i(t)$','interpreter','latex')
legend('Prey $[h=10^{-2}]$','Predator $[h=10^{-2}]$',...
    'Prey $[h=10^{-4}]$','Predator $[h=10^{-4}]$','interpreter','latex')
title({'Forward Euler with $h=10^{-2},h=10^{-4}$',...
    'Function Plot, $y_1(0)=100,y_2(0)=10$'},'interpreter','latex')
figure
plot(y1{1},y2{1},y1{2},y2{2})
xlabel('$y_1(t)$','interpreter','latex')
ylabel('$y_2(t)$','interpreter','latex')
legend('$y_1$ vs. $y_2$ $[h=10^{-2}]$','$y_1$ vs. $y_2$ $[h=10^{-4}]$',...
    'interpreter','latex')
title({'Forward Euler with $h=10^{-2},h=10^{-4}$',...
    'Phase Portrait, $y_1(0)=100,y_2(0)=10$'},'interpreter','latex')    
%% Problem 5
earthRadius = 6.378*10^6; % Radius of earth
moonRadius = 1.7374*10^6; % Radius of moon
d = 4.669*10^6; % [m]
D = 3.844*10^8; % [m]
tspan = [0 2.4*10^6];
x0 = 4.613*10^8; xdot0 = 0;
y0 = 0; ydot0 = -1074;
u0 = [x0,xdot0,y0,ydot0];
options = odeset('RelTol', 1e-6); 
[t,u] = ode45(@odefunc5,tspan,u0,options);
distance = zeros(length(t),1);
for j = 1:length(t)
    % getting euclidean distance from origin
    distance(j) = nthroot(((u(j,1)+d)^2 + u(j,3)^2),2); 
end
DistanceFromEarth = min(distance(:)) - earthRadius;
h = zeros(length(t)-1,1);
flag = 1;
for i = 1:length(t)-1
    h(i) = t(i+1)-t(i);
end
figure
plot(t(2:end),h)
xlabel('Time','interpreter','latex')
ylabel('Step Size $h$','interpreter','latex')
title({'Time Step Size, $h$, Evolution','Error Tolerance,$tol=10^{-6}$'}...
    ,'interpreter','latex')
figure
hold on
earth = nsidedpoly(1000, 'Center', [0 0], 'Radius', earthRadius);
moon = nsidedpoly(1000, 'Center', [D 0], 'Radius', moonRadius);
plot(earth, 'FaceColor', 'blue')
hold on
text(-35*d,0,'Earth $\rightarrow$','interpreter','latex')
plot(moon, 'FaceColor', 'cyan')
text(D-30*d,0,'Moon $\rightarrow$','interpreter','latex')
plot(u(:,1),u(:,3));
xlabel('X','interpreter','latex')
ylabel('Y','interpreter','latex')
legend('Earth','Moon','Sattelite Trajectory','interpreter','latex')
title({'Sattelite Trajectory','Computed using $\texttt{ode45}$',...
    'Error Tolerance, $tol=10^{-6}$'},'interpreter','latex')

%% Problem 5 - INTERNET FOUND SOLUTION [FOUND AFTER WRITING MY OWN](similar to mine, but not quite the same)
clear all

clc

% function cp09_08 %   three-body orbit problem
d = 4.669e6;
y0=[4.613e8; 0;  0; -1074];
time = zeros(6);
for k = 1:6    
    tol = 10^(-k);  
    options = odeset('RelTol', tol,  'AbsTol', tol*1e-3);
    tic;
    [t,y] = ode45(@apollo,[0,2.4e6], y0,  options);
    time(k) = toc;
end
figure;
    subplot(2,1,1);
    plot(t, y(:,1));
    xlabel('t');
    ylabel('x');
    title('Computer Problem 9.8  -- Solution of  ODE  for Three-Body Problem');
    subplot(2,1,2);
    plot(t, y(:,3));
    xlabel('t');
    ylabel('y');
    figure;
    s = t(2:end)-t(1:end-1);
    plot(t(1:end-1),s);
    xlabel('t');
    ylabel('step size');
    title('Computer Problem 9.8  -- Step Size  for Three-Body Problem');
        figure;
        plot(y(:,1), y(:,3));
        axis  equal;
        xlabel('x');
        ylabel('y');
        hold  on;
        title('Computer Problem 9.8  -- Computed  Orbit for Three-Body Problem');
            plot(-d,0,'o');
            plot(3.844e8-d,0,'o')
            text(-3d7,-5e7,'earth');
            text(3.5e8,-5e7,'moon')
            disp('Minimum   distance from  spacecraft earths surface:');
            min_dist = sqrt(min((y(:,1)+d).^2+y(:,3).^2))-6.378e6
            disp('Time to  compute solution for each tolerance:');
                disp('	tol 	time');
                for k = 1:6    
                    tol = 10^(-k);
                    fprintf('%8.1e  %12.7f\n', tol, time(k));
                end;

                disp(' ');
                   function [yprime] = apollo(t,y);
%                    y(1) = x,  y(2) = x', y(3) = y,  y(4) = y',G = 6.67259e-11;
%                    x = y(1); x' = y(2); y = y(3); y' = y(4); 
                   G = 6.67259e-11;
                   M = 5.974e24;
                   m = 7.348e22;
                   must = M/(m+M);
                   mu = m/(m+M);
                   D = 3.844e8;
                   d = 4.669e6;
                   omega = 2.661e-6;
                   r1 = sqrt((y(1)+d)^2+y(3)^2);
                   r2 = sqrt((D-d-y(1))^2+y(3)^2);
                   yprime = [y(2); -G*(M*(y(1)+mu*D)/r1^3+m*(y(1)-must*D)/r2^3)+omega^2*y(1)+...
                       2*omega*y(4); y(4);  -G*(M*y(3)/r1^3+m*y(3)/r2^3)+omega^2*y(3)-2*omega*y(2)];
                   end
% end
%% Functions
function [y_next] = ForwardEuler(func,h,y_previous)
y_next = y_previous + h*func;
end

function [y_next] = myfunc2(y_previous)
g = 21.937;
y_next = g*(1-(y_previous/80)^(3/2));
end

function dydt = predprey4(t,y)
alpha1 = 1; alpha2 = 0.5;
beta1 = 0.1; beta2 = 0.02;
dydt = zeros(2,1);
dydt(1) = y(1)*(alpha1-beta1*y(2));
dydt(2) = y(2)*(-alpha2 + beta2*y(1));
end

function dudt = odefunc5(t,u)
x = u(1);
xdot = u(2);
y = u(3);
ydot = u(4);
% Define Variables
G = 6.67259*10^(-11); % [m^3/(kg*s^2)]
M = 5.974*10^24; % [kg]
m = 7.348*10^22; % [kg]
mustar = M/(m+M);
mu = m/(m+M);
D = 3.844*10^8; % [m]
d = 4.669*10^6; % [m]
r1 = ((x+d)^2 + y^2)^(1/2);
r2 = ((D-d-x)^2 + y^2)^(1/2);
omega = 2.661*10^(-6); % [1/s]==[Hz]
xddot = -G*((M*(x+mu*D)/(r1^3))+m*(x-mustar*D)/(r2^3))...
    + (omega^2)*x + 2*omega*ydot;
yddot = -G*((M*y)/(r1^3) + ((m*y)/r2^3)) + (omega^2)*y...
    - 2*omega*xdot;
dudt = [xdot;xddot;ydot;yddot];
end
%%

                

                        
                        
                        
                        
                        
                        