%% start
% Julian days conversion
t1 = date2mjd2000([2016, 3, 14, 12, 0, 0]); %earth time
t2 = date2mjd2000([2016, 10, 15, 12, 0, 0]); %mars time

% Positions
[r1,v1] = EphSS_car(3,t1);
[r2,v2] = EphSS_car(4,t2);

muSun = getAstroConstants('Sun','mu');

%% langrangian coefficients function
function rf = FGKepler_trA(r0,v0,th, mu)

    h = cross(r0,v0);     % angolar momentum
    p = norm(h)^2 / mu; % semilatus rectum
    sigma0 = (dot(r0,v0)/sqrt(mu));
    R0 = norm(r0);
    Rf = (p*R0) / (R0 + (p-R0)*cos(th) - sqrt(p)*sigma0*sin(th));
    F = 1 - (Rf/p)*(1-cos(th));
    G = ((Rf*norm(r0)) / sqrt(mu*p)) * sin(th);
    rf = F*r0 + G*v0;

end


%% lambert arc with velocity output
function v1 = LMinETransfer(r1,r2,tm,mu)
    
    % position norms
    R1 = norm(r1);
    R2 = norm(r2);

    % semi-minor axis c
    c = norm(r2-r1);

    % semilatus rectum ro
    costh = (dot(r1,r2))/(R1 * R2);
    pmin = ((R1*R2)/c) *(1-costh);
    
    F = 1 - (R2/pmin)*(1-costh);

    sinth = tm*sqrt(1-costh^2);
    G = (R2*R1)/sqrt(mu*pmin) * sinth;

    v1 = 1/G * (r2-F*r1);
end

%% transfer plot
tm = 1;
vsc = LMinETransfer(r1,r2,tm,muSun);

% Plot transfer
nPoints=100; % number of points to plot for each orbit.
figure
hold on
title('Transfer Plot')
grid on
axis equal

dtheta =linspace (0,2*pi,nPoints); % theta vector

rSc_motion= zeros(nPoints, 3); % Initialize 

for iT=1:nPoints
    rSc_motion (iT, :)=FGKepler_trA(r1, vsc, dtheta(iT), muSun) ;
end
h_transfer_discrete = plot3(rSc_motion (:,1), rSc_motion (:,2), rSc_motion (:,3), Color= 'Red', Marker='*');


%% compare transfer orbit with ODE solver
%Equation of motion "r+mu*r/R^3 = 0
F=@ (t,x) [x(4); %dx/dt=Vx
    x(5); %dy/dt=Vy
    x(6); %dz/dt=Vz
    -muSun*x(1) / (sqrt (x(1)^2 +x(2)^2 +x(3)^2)^3) ; %dVx/dt
    -muSun*x(2) / (sqrt (x(1)^2 +x(2)^2 +x(3)^2)^3) ; %dVy/dt
    -muSun*x(3) / (sqrt (x(1)^2 +x(2)^2 +x(3)^2)^3)]; %dVz/dt

day = 86400; %seconds in a day
tspan = [0,(t2-t1)*day];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
y0 = [r1 v1];
[t,y] = ode45(F,tspan,y0,options); % ode solver

h_transfer_ode = plot3(rSc_motion(:, 1), rSc_motion(:, 2), rSc_motion(:, 3), ...
      'Color', [0, 1, 0, 0.5], 'LineWidth', 3); % RGBA color (green with 50% opacity)

%% orbit plot

% sun surface
%plot_planet(0, 0, 0, 5*1000*6378, 'ATATD-Toolbox/TButils/textures/Sun.png')
plot_planet(0, 0, 0, 5*1000*6378, 'dani.jpg')


hold on
% Plot Earth

% earth surface
plot_planet(r1(1), r1(2), r1(3), 3*1000*6378, 'ATATD-Toolbox/TButils/textures/Earth.jpg')

rEarth_motion= zeros(nPoints, 3); % Initialize 

for iT=1:nPoints
    rEarth_motion (iT, :)=FGKepler_trA(r1, v1, dtheta(iT), muSun) ;
end
h_earth_orbit = plot3(rEarth_motion (:,1), rEarth_motion (:,2), rEarth_motion (:,3), LineWidth=3);


% Plot mars

% mars surface
plot_planet(r2(1), r2(2), r2(3), 3*1000*3390, 'ATATD-Toolbox/TButils/textures/Mars.jpg')

rMars_motion= zeros(nPoints, 3); % Initialize 

for iT=1:nPoints
    rMars_motion (iT, :)=FGKepler_trA(r2, v2, dtheta(iT), muSun) ;
end
h_mars_orbit = plot3(rMars_motion (:,1), rMars_motion (:,2), rMars_motion (:,3),  LineWidth=3 );

% Create legend using plot handles
legend([h_transfer_discrete, h_transfer_ode, h_earth_orbit, h_mars_orbit], ...
       {'Transfer Orbit (Discrete Solution)', ...
        'Transfer Orbit (ODE45 Solution)', ...
        'Earth Orbit', ...
        'Mars Orbit'}, ...
       'Location', 'best');
