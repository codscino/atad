
%% start
% Julian days conversion
t1 = date2mjd2000([2016, 03, 14, 0, 0, 0]); %earth time
t2 = date2mjd2000([2016, 10, 15, 0, 0, 0]); %mars time

% ex. 6 extended dates (it crashes)
%t1 = date2mjd2000([2016, 09, 1, 12, 0, 0]); %earth time
%t2 = date2mjd2000([2016, 11, 19, 12, 0, 0]); %mars time

% Positions
[r1,v1] = EphSS_car(3,t1);
[r2,v2] = EphSS_car(4,t2);

muSun = getAstroConstants('Sun','mu');

% initial velocity of the spacecraft
tm = 1;
vsc = LMinETransfer(r1,r2,tm,muSun);



%% first iteration
% find Smat
dT = 215*86400;
Smat = STM_Lambert(r1, vsc, dT, muSun);
[rSc_final, vSc_final] = FGKepler_dt2(r1, vsc, dT, muSun) ;
dr_t2 = r2 - rSc_final;


%% while loop
Error = 10;
tol = 10^(-3); %1 metre tolerance
numIter = 0;

while Error > tol
    numIter = numIter + 1;
    Smat = STM_Lambert(r1, vsc, dT, muSun);
    dv_t1 = inv(Smat) * dr_t2';
    vsc = vsc + dv_t1';
    [rSc_final, vSc_final] = FGKepler_dt2(r1, vsc, dT, muSun) ;
    dr_t2 = r2 - rSc_final;
    Error = norm(dr_t2);
end

%% compute total deltav
vSc_initial = vsc; 

deltav1 = norm(v1-vSc_initial);
deltav2 = norm(v2-vSc_final);
deltavTot = deltav1 + deltav2