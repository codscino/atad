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



%% first iteration
% find Smat
dT = 215*86400;
Smat = STM_Lambert(r1, vsc, dT, muSun);
rSc_final = FGKepler_dt(r1, vsc, dT, muSun) ;
dr_t2 = r2 - rSc_final;


%% while loop
niteration = 0;
Error = 10;

while Error > 10^(-3) %1 metre
    niteration = niteration + 1;
    Smat = STM_Lambert(r1, vsc, dT, muSun);
    dv_t1 = inv(Smat) * dr_t2';
    vsc = vsc + dv_t1';
    rSc_final = FGKepler_dt(r1, vsc, dT, muSun) ;
    dr_t2 = r2 - rSc_final;
    Error = norm(dr_t2);
end
