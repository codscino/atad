%% start
% Julian days conversion
t1 = date2mjd2000([2016, 3, 14, 12, 0, 0]); %earth time
t2 = date2mjd2000([2016, 10, 15, 12, 0, 0]); %mars time

% Positions
[r1,v1] = EphSS_car(3,t1);
[r2,v2] = EphSS_car(4,t2);

muSun = getAstroConstants('Sun','mu');

%% lambert algorithm
function [amin, emin, dtmin] = MinETransfer(r1,r2,tm,mu)
    
    % position norms
    R1 = norm(r1);
    R2 = norm(r2);

    % semi-minor axis c
    c = norm(r2-r1);

    % semi-major axis
    amin = 1/4 * (R1+R2+c);

    % semilatus rectum ro
    costh = (dot(r1,r2))/(R1 * R2);
    pmin = ((R1*R2)/c) *(1-costh);
    
    % eccentricty
    emin = sqrt(1-(pmin/amin));

    % beta angle
    beta = 2*asin(sqrt((2*amin - c) / (2*amin)));
    dtmin = sqrt(amin^3/mu)*(pi-tm*(beta-sin(beta)));
end

%% final results
tm = 1;
[amin, emin, dtmin] = MinETransfer(r1,r2,tm,muSun);
days = dtmin/86400;