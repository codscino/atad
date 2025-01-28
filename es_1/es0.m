%% dani excercise
clc;
clear;

% 2 body problem diff equation
% r dot dot + mu_earth/r^3 * r = 0

y0 = [-18676 6246 12474 0.551 -1.946 -3.886];
mu = 3.986e5;

% THE EARTH AROUND THE SUN
% y0 = [0 149597871 0 29.784 0 0];
% mu = 1.327e11;

day = 86400; %seconds in a day
tspan = [0,1*day];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

% convert to a I order ODE with v and vdot
dydt = @(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)];  % [v,vdot]
[t,y] = ode45(dydt,tspan,y0,options); % ode solver

%plot the ode solver result
figure;
plot3(y(:,1), y(:,2), y(:,3));
axis equal;
grid on;
xlabel('km');
ylabel('km');
zlabel('km');
title('Codina','FontSize', 20);

