% Author: Nicholas Piercy
% Date : 12/12/2021
% Numerical Methods, Homework #8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% Problem 1
x0 = 1; x2 = 5; % boundary conditions
dxdt1 = 1; dxdt2 = 1.6; % initial guesses
InitCond1 = [x0,dxdt1]; InitCond2 = [x0,dxdt2];
stepSize = 0.25;
tspan = [0:stepSize:2];
[t1,u1] = ode45(@odeP1,tspan,InitCond1);
[t2,u2] = ode45(@odeP1,tspan,InitCond2);
alpha = (x2-u2(end,1))/(u1(end,1)-u2(end,1));
w = alpha*u1(:,1)+(1-alpha)*u2(:,1);
plot(t1,u1(:,1))
hold on
plot(t2,u2(:,1))
plot(tspan,w)
title({'BVP Solution using the Shooting Method',['Step Size t = ',num2str(stepSize)],...
    '$\frac{d^2x(t)}{dt^2} + t \frac{dx(t)}{dt} - 3x(t) = 3t$','$ x(0) = 1, x(2) = 5$'}...
    ,'Fontsize',15,'interpreter','latex')
legend(['$\frac{dx(t)}{dt}_1 = $',num2str(dxdt1)],['$\frac{dx(t)}{dt}_2 = $'...
    ,num2str(dxdt2)],'$w(t)$','location','best','Fontsize',20,'interpreter','latex')
xlabel('time (t)','Fontsize',20,'interpreter','latex')
ylabel('x(t)','Fontsize',20,'interpreter','latex')

%% Problem 2
numTimeSteps = 9;
interval = [0,4];
h = length(interval)/numTimeSteps;
x = linspace(interval(1),interval(2),numTimeSteps); % discretize function space
LHS = zeros(numTimeSteps-2,numTimeSteps-2);
RHS = x';
for i = 1:numTimeSteps-2 % only interior points
    LHS(i,i) = -(4+6*h^2)/(2*h^2);
    if i < numTimeSteps-2; LHS(i,i+1) = (2-h)/(2*h^2); end
    if i > 1; LHS(i,i-1) = (2-h)/(2*h^2); end
end
disp(LHS)
% account for boundary terms
LeftBoundary = 0;
RightBoundary = 5;
RHS(2) = RHS(2)-LeftBoundary*(2-h)/(2*h^2); 
RHS(end-1) = RHS(end-1)-RightBoundary*(2-h)/(2*h^2);

[Thom_U,Thom_L] = ThomasAlgorithm(LHS);
soln = BackwardSub(Thom_U,RHS(2:numTimeSteps-1));
soln = [0;soln;5];
plot(x,soln)
hold on
for i = 1:numTimeSteps
    exact(i) = 6.05722*10^(-8)*(2.71828^(14.4222 - 2.30278*x(i))+...
    3*x(i) - 58*2.71828^(-2.30278*(x(i) - 4)) - 2.71828^(1.30278*x(i))...
        + 58*2.71828^(1.30278*x(i) + 9.2111) - 1.83436*10^6*(3*x(i) + 1) + 1);
end
plot(x,exact)
title({'$\frac{d^2u}{dx^2} + \frac{du}{dx} - 3u = x$',...
    ['Number of discretization points = ',num2str(numTimeSteps)]}...
    ,'Fontsize',20,'interpreter','latex')
xlabel('x','Fontsize',20,'interpreter','latex')
ylabel('y','Fontsize',20,'interpreter','latex')
legend('FDM','exact','Fontsize',20,'interpreter','latex')

%% Problem 3 [Computer Proplem 10.1a]
x0 = 0; x2 = 1; % boundary conditions
% Discretization
stepSize = 0.05;
tspan = [0 1];
numTimeSteps = ((tspan(2)-tspan(1))/stepSize);
t = linspace(tspan(1),tspan(2),numTimeSteps);
% Stopping criterion
maxIterations = 50;
RelErrorTol = 1e-5;
dudt(1) = 1; %initial guess 
InitCond1 = [x0,dudt(1)];
u{1} = ode4(@odeP3,tspan(1),stepSize,tspan(2),InitCond1);
plot(t,u{1}(:,1))
hold on
%Guess 2
dudt(2) = 2; % initial guess
InitCond2 = [x0,dudt(2)];
u{2} = ode4(@odeP3,tspan(1),stepSize,tspan(2),InitCond2);
plot(t,u{2}(:,1))
%Linear Combination to find next guess
alpha1 = (x2-u{2}(end,1))/(u{1}(end,1)-u{2}(end,1));
dudt(3) = alpha1*dudt(1)+(1-alpha1)*dudt(2);  
for i = 3:maxIterations
    InitCond = [x0, dudt(i)];
    u{i} = ode4(@odeP3,tspan(1),stepSize,tspan(2),InitCond);
    plot(t,u{i}(:,1))
    alpha = (x2-u{i}(end,1))/(u{i-1}(end,1)-u{i}(end,1));
    error = abs(u{i}(end,1) - x2);
    RelativeError = abs(u{i}(end,1)-u{i-1}(end,1));
    if RelativeError <= RelErrorTol
        break
    elseif error <= RelErrorTol
        break
    end
    dudt(i+1) = alpha*dudt(i-1) + (1-alpha)*dudt(i); 
end
NumIterationToConverge = i;
w = alpha*u{NumIterationToConverge}(:,1)+(1-alpha)*u{NumIterationToConverge-1}(:,1);
plot(t,w,'linewidth',2)
Legend = cell(NumIterationToConverge+1,1);
for k=1:NumIterationToConverge
    Legend{k}=strcat('Iteration', num2str(k));
end
Legend{end} = ('Final Solution');
legend(Legend,'fontsize',20,'location','best','interpreter','latex')
title({'BVP Solution for Non-Linear ODE using the Shooting Method',...
    ['$u"(t) = 10u^3(t) + 3u(t) + t^2$',';      $ u(0) = 0, u(1) = 1$'],...
    ['Initial Guesses $z_1=$',num2str(dudt(1)),',$z_2=$',num2str(dudt(2)),...
    [',   with step Size t = ',num2str(stepSize)],]},'Fontsize',15,'interpreter','latex')
xlabel('time,t','Fontsize',20,'interpreter','latex')
ylabel('u(t)','Fontsize',20,'interpreter','latex')
message1 = ("Number of iterations to converge:");
message2 = ('Final Slope Value determined by Shooting Method');
disp(message1)
disp(num2str(NumIterationToConverge))
disp(message2)
disp(num2str(dudt(end)))

%% Problem 3 [Computer Proplem 10.1b]
clear all
close all
clc

numSteps = [3 7 15];
tspan = [0 1];
y0 = [0 1]; % initial condition
for i = 1:length(numSteps)
    n = numSteps(i);
    h = 1/(n-1);
    t = linspace(tspan(1),tspan(2),n);
    
    u = zeros(n,1);
    u(1) = y0(1);
    u(n) = y0(2);
    
%         lhs = zeros(n-2,n-2);
%         rhs = zeros(n-2,1);
    rhs_func = @(t,u) 10*u.^3 + 3*u;
    func = rhs_func(t,u);
    for j = 1:n-2
        if j ==1
            rhs(j) = 3*t(j+1) - u(1)/h^2;
%             lhs(j,j+1) = (-2*u(j)+u(j-1))/h^2 - func; 
            lhs(j,j) = -2/h^2 - func;
            lhs(j,j+1) = 1/h^2;
        elseif j == n-2 
            rhs(j) = 3*tval - u(2)/h^2;
%             lhs(j,j) = (u(j+1)-2*u(j))/h^2 - func;
            lhs(j,j-1) = -2/h^2 - func;
            lhs(j,j) = 1/h^2;
        else
            rhs(j) = 3*t(j+1);
%             lhs(j,j-1) = (u(j+1)-2*u(j)+u(j-1))/h^2 - func;
            lhs(j,j-1) = 1/h^2;
            lhs(j,j) = -2/h^2 - func;
            lhs(j,j+1) = 1/h^2;
        end
    end
    [Thom_U,Thom_L] = ThomasAlgorithm(lhs);
    soln = BackwardSub(Thom_U,rhs);
    soln = [y0(1);soln;y0(2)];
    plot(t,soln)
    hold on
    end

    
    
    
    
    



%%

syms u
LeftBoundary = 0;
RightBoundary = 1;
numTimeSteps = [10,3,7,15];
for k = 1:length(numTimeSteps)
    n = numTimeSteps(k);

tspan = [0,1];
h = (tspan(2)-tspan(1))/(n-1);
t = linspace(tspan(1),tspan(2),n-1); % discretize function space
LHS = zeros(n-2,n-2);
y0 = [0, 1];

% RHS = @(t,y) 10*y^3+3*y+t^2;
soln = finiteDiff(y,t,y0,n);
soln{k} = [y0(1); soln; y0(2)];
end
function LHS = finiteDiff(y)

for i = 2:n-1 % only interior points
    RHS = @(t,y) 10*y(i)^3+3*y(i)+t(i)^2;
    LHS(i,i) = -2*y(i)/(h^2)-RHS(i);
    if i < n-2
        LHS(i,i+1) = y(i+1)/(h^2)-RHS(i); 
    end
    if i > 1
        LHS(i,i-1) = y(i-1)/(h^2)-RHS(i); 
    end
%     y(i) =  RHS
end
% disp(LHS)
% account for boundary terms
% for i = 1:n-2
%     RHS(i) = @(t,y) 10*y(i)^3+3*y(i)+t(i)^2;
% end
% rhs = @(u,x) (10*u.^3 + 3.*u + x.^2);
% RHS = sym(rhs,[n-2,1]);
% RHS(2) = RHS(2)-LeftBoundary*1/(h^2); 
% RHS(end-1) = RHS(end-1)-RightBoundary*1/(h^2);

% [Thom_U,Thom_L] = ThomasAlgorithm(LHS);
% soln = BackwardSub(Thom_U,RHS(2:n-1));
% soln = [LeftBoundary;soln;RightBoundary];
% plot(x,soln)
% hold on
% for i = 1:numTimeSteps
%     exact(i) = 6.05722*10^(-8)*(2.71828^(14.4222 - 2.30278*x(i))+...
%     3*x(i) - 58*2.71828^(-2.30278*(x(i) - 4)) - 2.71828^(1.30278*x(i))...
%         + 58*2.71828^(1.30278*x(i) + 9.2111) - 1.83436*10^6*(3*x(i) + 1) + 1);
% end
% plot(x,exact)
end
% title({'$\frac{d^2u}{dx^2} + \frac{du}{dx} - 3u = x$',...
%     ['Number of discretization points = ',num2str(n)]}...
%     ,'Fontsize',20,'interpreter','latex')
% xlabel('x','Fontsize',20,'interpreter','latex')
% ylabel('y','Fontsize',20,'interpreter','latex')
% legend('FDM','exact','Fontsize',20,'interpreter','latex')
%%
% % % % % % %% Problem 3 - chegg
% % % % % % clear all
% % % % % % close all
% % % % % % clc
% % % % % % % function cp10_01 % two-point boundary value problem
% % % % % % global a b ua ub rhs ivp;
% % % % % % a = 0;
% % % % % % b = 1;
% % % % % % tspan = [a,b];
% % % % % % ua = 0;
% % % % % % ub = 1;
% rhs = inline('10*u.^3+3*u+t.^2','t','u');
% rhs = @(t,u) 10*u.^3+3*u+t.^2';
% % % % % % % ivp = inline('[y(2);10*y(1)^3+3*y(1)+t^2]','t','y');
% % % % % % ivp = @(t,y) [y(2);10*y(1)^3+3*y(1)+t^2];
% % % % % % disp('(a) Initial slope determined by shooting method:');
% % % % % % slope = fzero(@shooting_fun,0.5);
% % % % % % [t, y] = ode45(ivp,[a,b],[ua;slope]);
% % % % % % plot(t, y(:,1),'r-');
% % % % % % xlabel('t');
% % % % % % ylabel('u');
% % % % % % axis([0 1 0 1.2]);
% % % % % % title('Computer Problem 10.1(a) -- Shooting Method for Two-Point BVP');
% n = [1 3 7 15];
% figure;
% options = optimset('MaxFunEvals', 5000, 'Display', 'off');
% for i = 1:4
%     t = linspace(a,b,n(i)+2);
%     y = [ua; fsolve(@fd_fun,t(2:n(i)+1),options)'; ub];
%     plot(t,y);
%     hold on;
% end
% % % % % % xlabel('t');
% % % % % % ylabel('u');
% % % % % % axis([0 1 0 1.2]);
% % % % % % title('Computer Problem 10.1(b) -- Finite Difference Method for Two-Point BVP');
% % % % % % % % % % % % n = [3 4 5 6];
% % % % % % % % % % % % figure;
% % % % % % % % % % % % disp('(c) Coefficients of polynomials determined by collocation:');
% % % % % % % % % % % % for i = 1:4
% % % % % % % % % % % %     x = fsolve(@colloc_fun,ones(n(i),1),options);
% % % % % % % % % % % %     t = linspace(a,b,50);
% % % % % % % % % % % %     plot(t,polyval(x,t));
% % % % % % % % % % % %     hold on;
% % % % % % % % % % % % end
% % % % % % % % % % % % xlabel('t');
% % % % % % % % % % % % ylabel('u');
% % % % % % % % % % % % axis([0 1 0 1.2]);
% % % % % % % % % % % % title('Computer Problem 10.1(c) -- Collocation Method for Two-Point BVP');
% % % % % % %% Functions
% % % % % % function [g] = shooting_fun(slope);
% % % % % %     global a b ua ub ivp;
% % % % % %     [t, y] = ode45(ivp,[a,b],[ua;slope]);
% % % % % %     plot(t, y(:,1));
% % % % % %     hold on;
% % % % % %     g = ub-y(end,1);
% % % % % % end
% function [z] = fd_fun(y);
%     global a b ua ub rhs;
%     n = length(y);
%     h = 1/(n+1);
%     t = linspace(h,1-h,n);
%     f = rhs(t,y);
%     if n>1
%         z = [(y(2)-2*y(1)+ua)/h^2-f(1); ((y(3:n)-2*y(2:n-1)+y(1:n-2))/h^2-f(2:n-1))';(ub-2*y(n)+y(n-1))/h^2-f(n)];
%     else
%         z = (ub-2*y+ua)/h^2-f;
%     end
% end
% % % % % % function [z] = colloc_fun(x);
% % % % % % % collocation method
% % % % % %     global a b ua ub rhs;
% % % % % %     n = length(x);
% % % % % %     h = 1/(n-1);
% % % % % %     t = linspace(a,b,n)';
% % % % % %     d = (n-1:-1:1)'.*x(1:n-1);
% % % % % %     d2 = (n-2:-1:1)'.*(d(1:n-2));
% % % % % %     z = [polyval(x,a)-ua; polyval(d2,t(2:n-1))-rhs(t(2:n-1), polyval(x,t(2:n-1))); polyval(x,b)-ub];
% % % % % % %     z = [feval(x,a)-ua; feval(d2,t(2:n-1))-rhs(t(2:n-1), feval(x,t(2:n-1))); feval(x,b)-ub];
% % % % % % end
% % % % % % 
% % % % % % 
% % % % % %                    
% % % % % %                 
% % % % % %             


%% Functions
function dxdt = odeP1(t,u)
x = u(1);
xdot = u(2);
xddot = -t*xdot + 3*x -3*t;
dxdt = [xdot;xddot];
end

function dudt = odeP3(t,u)
x = u(1);
xdot = u(2);
xddot = 10*x^3 + 3*x + t^2;
dudt = [xdot;xddot];
end

function yout = ode4(F,t0,h,tfinal,y0)
y = y0;
yout = y0;
for t = t0+h : h : tfinal-h
     s1 = F(t,y);
     s2 = F(t+h/2, y+h*s1/2);
     s3 = F(t+h/2, y+h*s2/2);
     s4 = F(t+h, y+h*s3);
     % solution for first and second derivatives
     yval = y(1) + h*(s1(1) + 2*s2(1) + 2*s3(1) + s4(1))/6; 
     yprime = y(2) + h*(s1(2) + 2*s2(2) + 2*s3(2) + s4(2))/6;
     % new initial vector for next time step
     y = [yval,yprime];
     yout = [yout; y];
end   
end



