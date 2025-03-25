% VortexPanelMethod.m
% This script uses the Vortex Panel Method to calculate the pressure coefficient of NACA 4412 airfoil

% Load airfoil coordinates
run('AirfoilNACA4412.m');

%% Combine coordinates of the upper and lower surfaces
x_airfoil = [xu(1:end-1), flip(xl(1:end))]; % Combine upper and lower surfaces
y_airfoil = [yu(1:end-1), flip(yl(1:end))];

%% Panel method setup
N = length(x_airfoil) -1; % Number of panels
panels = struct('x',{},'y',{},'xp',{},'yp',{},'theta',{},'gamma',{});
% To make sure the airfoil is closed
% figure;
% plot(x_airfoil,y_airfoil,'b','LineWidth',2);
% axis equal;
% hold on;

% Define panels
for i = 1:N
    panels(i).x = [x_airfoil(i), x_airfoil(i+1)];
    panels(i).y = [y_airfoil(i), y_airfoil(i+1)];
    panels(i).xp = (panels(i).x(1) + panels(i).x(2))/2; % Panel center x-coordinate
    panels(i).yp = (panels(i).y(1) + panels(i).y(2))/2; % Panel center y-coordinate
    dx = panels(i).x(2) - panels(i).x(1);
    dy = panels(i).y(2) - panels(i).y(1);
    panels(i).theta = atan2(dy,dx); % Panel angle
    panels(i).gamma = 0;         % Initialize vortex strength to zero
end

%% Free stream vilecity
V_inf = 100;
alpha = 0;
alpha_rad = deg2rad(alpha);

%% Calculate free stream velocity for each panel
b = zeros(N,1);
for i = 1:N
    % Calculate normal component of free stream velocity at panel center
    b(i) = V_inf * cos(alpha_rad - panels(i).theta);
end

%% Influence coefficient matrix
A = zeros(N,N);
for i = 1:N
    for j = 1:N
        if i ~= j
            % Calculate induced velocity at panel i due to vortex on panel j
            [u,v] = induced_velocity(panels(j).xp,panels(j).yp,panels(j).gamma,panels(i).xp,panels(i).yp);
            u_normal = u*cos(panels(i).theta) + v*sin(panels(i).theta);
            A(i,j) = u_normal;
        else
            % Diagonal term for Kutta condition
            A(i,j) = 0.5;
        end
    end
end

%% Kutta condition
A(end,:) = 0;
A(end,end) = 1;
A(end,1) = 1;
b(end) = 0;

%% Solve for vortex strength
gamma = A \ b;

%% Calculate velocity diistribution
%   V = zeros(N,1);
%   for i = 1:N
%        V(i) = 0;
%       for j = 1:N
%           [u,v] = induced_velocity(panels(j).xp,panels(j).yp,gamma(j),panels(i).xp,panels(i).yp);
%           V(i) = V(i) + u*cos(panels(i).theta) + v*sin(panels(i).theta);
%      end
%      V(i) = V(i) + V_inf*cos(alpha_rad + panels(i).theta);
%      % Add free stream velocity
%  end 
V = -1 .* gamma;

%% Calculate pressure coefficient
Cp = 1 - (gamma/V_inf).^2;

%% Plot results
figure;
plot([panels.xp], -Cp, 'o-');
xlabel('x/c');
ylabel('Cp');
title('Pressure coefficient distribution');

%% Function for induced velocity
function [u,v] = induced_velocity(xv,yv,gamma,x,y)
    % Calculate insuced velocity at point (x,y) due to vortex at (xv,yv)
    dx = x - xv;
    dy = y - yv;
    r = sqrt(dx^2 + dy^2);
    u = -gamma/(2*pi*(r^2))*dy;
    v = gamma/(2*pi*(r^2))*dx;
end
