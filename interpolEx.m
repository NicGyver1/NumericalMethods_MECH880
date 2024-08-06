% Runge example of problems with the interpolating 
% polynomial (high oscillations) when using
% more and more points to interpolate through,
% and how to fix it using piecewise polynomial interpolation
% (cubic-splines)

clf

% input data
n1=6; % number of interpolation points 
n2=21; % a larger number of interpolation points

% create two vectors with the x-coordinates of n1, n2 equally spaced values
% inside the [-5,5] interval
x1=linspace(-5,5,n1);
x2=linspace(-5,5,n2);

% the corresponding y-values for the Runge function f(x)=1/(1+x^2); see the
% function called "givenf" after the script 
y1=givenf(x1);
y2=givenf(x2);

% number of points where we will plot the various polynomials, functions.
nplot=200;

x=linspace(-5,5,nplot);

% the y-values for the Runge function at the plotting points
ygiven=givenf(x);

% the y-values for the interpolating polynomial (using Lagrange's form, see
% function at the end)
ylag1 = zeros(1,nplot);
ylag2 = zeros(1,nplot);
for i=1:nplot
    ylag1(i)=lagrange(n1,x(i)); % for the polynomial through n1 points
    ylag2(i)=lagrange(n2,x(i)); % for the polynomial through n2 points
end

nw=3; %linewidth
figure(1) % plotting 
% create text to show in the Legend

str1 = [num2str(n1), ' equally-spaced data points'];
str2 = [num2str(n2), ' equally-spaced data points'];
first_Lagrange_pol = ['Lagrange interp. pol. of degree ', num2str(n1-1)];
first_cubic_spline = ['cubic spline for ', num2str(n1), ' data points'];
second_Lagrange_pol = ['Lagrange interp. pol. of degree ', num2str(n2-1)]; 
second_cubic_spline = ['cubic spline for ', num2str(n2), ' data points'];


plot(x,ygiven,'Linewidth', nw) % plotting the Runge function
axis([-5 5 -.5 2]) 
legend('The Runge function')
hold on
pause  % press any key to continue

plot(x1,y1, 'go', 'Linewidth', nw, 'markerSize', 14) % the interpolating data for the first set
legend('The Runge function', str1)
pause % press any key to continue

plot(x,ylag1,'g--', 'Linewidth', nw) % the first interpolating polynomial
legend('The Runge function', str1, first_Lagrange_pol)
pause  % press any key to continue

yy=spline(x1,givenf(x1),x); % compute the cubic spline through the first set of data points at all the points where we plot it
plot(x,yy, 'm-.','Linewidth', nw);
legend('The Runge function', str1, first_Lagrange_pol,...
    first_cubic_spline)
pause % press any key to continue

plot(x2,y2, 'ko', 'Linewidth', nw, 'markerSize', 11)  % the interpolating data for the second set
legend('The Runge function', str1, first_Lagrange_pol,...
    first_cubic_spline, str2)
pause  % press any key to continue

plot(x,ylag2,'k:','Linewidth', nw) % the second interpolating polynomial
legend('The Runge function', str1, first_Lagrange_pol,...
    first_cubic_spline, str2, second_Lagrange_pol)
pause

yy=spline(x2,givenf(x2),x); % compute the cubic spline through the second set of data points
plot(x,yy, 'r','Linewidth', nw);

legend('The Runge function', str1, first_Lagrange_pol,...
    first_cubic_spline, str2, second_Lagrange_pol,...
    second_cubic_spline)



hold off
% 
% % example for non-shape preserving property of the cubic spline
% figure(2)
% hold on
% 
% x=linspace(-5, 5, 11);
% y=[1, 1.2, 1.3, 1.73, 3, 6, 15, 25 , 27, 85,2];
% xx=linspace(-5,5,200);
% yy=spline(x,y,xx);
% zz=pchip(x,y,xx);
% plot(x,y,'bo',xx,yy,'b', xx,zz,'g','Linewidth',nw)
% 
% y=[10, 2, 3, 9, 3, 14, 9, 1 , 15, 5,2];
% yy=spline(x,y,xx);
% zz=pchip(x,y,xx);
% 
% plot(x,y,'ro',xx,yy, 'r', xx,zz, 'b')
% 
% hold off

function [y]=givenf(x)
% The Runge example; 1/(1+x^2) over [-5,5]

y=1./(1+x.^2);
end

function y=lagrange(ndata,x)
% Lagrange interpolating polynomial of degree (ndata-1)
% returns the value of the polynomial at a point x
% uses givenf 

% number of data points =ndata
xdata=linspace(-5,5,ndata);
ydata=givenf(xdata);

sum=0;
for i=1:ndata
    prod=1;
    for j=1:ndata
        if i~=j 
            prod=prod*(x-xdata(j))/(xdata(i)-xdata(j));
        end
    end
    sum=sum+prod*ydata(i);
end
y=sum;
end