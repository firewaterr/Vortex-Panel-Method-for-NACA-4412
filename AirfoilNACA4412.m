clc;
clear;
close all;


%% Input parameter
m = 0.04;   % maximun camber
p = 0.4;    % location of maximum camber
t = 0.12;   % thickness


%% Calculate coordinates
n = 100;    % number of points on each side of the chord
x = linspace(0,1,n+1);  % x-coordinates
yt = 5*t*(0.2969*sqrt(x)-0.126*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4);  %thickness distribution
yf = zeros(size(x));
dyfdx = zeros(size(x));
% Calculate camber line
if p ~= 0    % Cambered airfoil
for i = 1:n+1
    if x(i)<=p
    yf(i) = (m/p^2)*(2*p*x(i)-x(i)^2);
    dyfdx(i) = (2*m/p^2)*(p-x(i));
    else
    yf(i)=(m/(1-p)^2)*((1-2*p)+2*p*x(i)-x(i)^2);
    dyfdx(i)=(2*m/(1-p)^2)*(p-x(i));
    end
end
end

% Calculate angle of the caber line
theta = atan(dyfdx);

% Calculate coordinates of the upper and lower surfaces;
xu = x - yt.*sin(theta);
xl = x + yt.*sin(theta);
yu = yf + yt.*cos(theta);
yl = yf - yt.*cos(theta);


%% Plot airfoil
figure;
plot(xu,yu,'b',xl,yl,'r','LineWidth',2);
axis equal;
xlabel('x/c');
ylabel('y/c');
legend('Upper surface','Lower surface');
title('NACA 4412 airfoil');
