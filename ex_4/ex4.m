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

% find Smat
dT = 215*86400;
Smat = STM_Lambert(r1, vsc, dT, muSun);

%% follow diagram flow
rSc_final = FGKepler_dt(r1, vsc, dT, muSun) ;
dr_t2 = r2 - rSc_final;

% norm(dr_t2) is bigger than the tolerance so i apply corrections
dv_t1 = inv(Smat) * dr_t2'; % velocity correction (it is a column vector)
vsc_new = vsc + dv_t1';

% try to re-propagate
rSc_final_new = FGKepler_dt(r1, vsc_new, dT, muSun) ;
dr_t2_new = r2-rSc_final_new;
norm(dr_t2_new)

