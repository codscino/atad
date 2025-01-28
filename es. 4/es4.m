%% start
% Julian days conversion
t1 = date2mjd2000([2016, 3, 14, 12, 0, 0]); %earth time
t2 = date2mjd2000([2016, 10, 15, 12, 0, 0]); %mars time

% Positions
[r1,v1] = EphSS_car(3,t1);
[r2,v2] = EphSS_car(4,t2);

muSun = getAstroConstants('Sun','mu');

% velocity of the spacecraft
tm = 1;
vsc = LMinETransfer(r1,r2,tm,muSun);


% STM(state transition matrix) lambert
function Smat = STM_Lambert(r0, v0, dt, mu)
    R0 = norm(r0);
    V0 = norm(v0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% this is basically FGKepler_dt function %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = mu / ((2*mu)/R0 - V0^2);
    n = sqrt(mu/a^3);
    dM = n*dt;
    sigma0 = dot(r0,v0) / sqrt(mu);

    % find dE with fzero
    fun = @(dE) -dM + dE - (1- (R0/a))*sin(dE) - (sigma0/sqrt(a))*(cos(dE)-1);
    x = fzero(fun,dM);
    dE = x;

    F = 1 - (a/R0)*(1-cos(dE));
    G = dt + sqrt(a^3/mu)*(sin(dE)-dE);

    rf = F*r0 + G*v0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% end of FGKepler_dt function %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Rf = norm(rf);
    dF = - sqrt(mu*a)/(Rf*R0) * sin(dE);
    dG = 1 - a / Rf *(1-cos(dE));
    vf = dF*r0 + dG*v0;

    C = a*sqrt(a^3/mu) * (3*sin(dE) - (2+cos(dE)) *dE) -a*dt*(1-cos(dE));
    deltar = rf-r0;
    deltav = vf-v0;
    
    % sensitivity matrix
    Smat = r0/mu * (1-F)*(deltar*v0' - deltav*r0') + ...
    (C/mu)*vf*v0' + G*eye(3);
end


% find Smat
dT = 215*86400;
Smat = STM_Lambert(r1, vsc, dT, muSun);

%% follow diagram
rSc_final = FGKepler_dt(r1, vsc, dT, muSun) ;
dr_t2 = rSc_final - r2;

% norm(dr_t2) is bigger than the tolerance so i apply corrections
dv_t1 = inv(Smat) * dr_t2; % velocity correction

rSc_initial = FGKepler_dt(r1, vsc, 0, muSun) ;
new_vsc = vsc + dv_t1;

% try to re-propagate
new_rSc_final = FGKepler_dt(r1, new_vsc, dT, muSun) ;
dr_t2 = new_rSc_final - r2;