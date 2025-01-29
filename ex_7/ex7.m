%% start
clear all
clc
% Julian days conversion
t1 = date2mjd2000([2016, 03, 14, 12, 0, 0]); % earth time
t2 = date2mjd2000([2016, 10, 15, 12, 0, 0]); % mars time

% Positions
[r1,v1] = EphSS_car(3,t1);
[r2,v2] = EphSS_car(4,t2);

muSun = getAstroConstants('Sun','mu');

% initial velocity of the spacecraft
tm = 1;
vsc = LMinETransfer(r1,r2,tm,muSun);
vSc_initial = vsc;


%% first iteration
% find Smat
dT = 215*86400; %target time of 215 days
Smat = STM_Lambert(r1, vsc, dT, muSun);
[rSc_final, vSc_final] = FGKepler_dt2(r1, vsc, dT, muSun) ;
dr_t2 = r2 - rSc_final;


%% while loop
Error = 10;
tol = 10^(-3); %1 metre tolerance
numIter = 0; % iterations while loops
N = 10; %total number of iterations for the for loop
[amin, emin, dtmin] = MinETransfer(r1,r2,tm,muSun); %dtmin from the min energy transfer

for i = 1:N %i index of the for loop
    lambda = i/N;
    dt = lambda*dT + (1-lambda)*dtmin;
    while (Error>tol) && (numIter<25)
        numIter = numIter + 1;
        [rSc_final, vSc_final] = FGKepler_dt2(r1, vsc, dT, muSun);
        Smat = STM_Lambert(r1, vSc_final, dT, muSun);
        dr_t2 = r2 - rSc_final;
        dv_t1 = inv(Smat) * dr_t2';
        vsc = vsc + dv_t1';
        Error = norm(dr_t2);
    end
end

%% compute total deltav
deltav1 = norm(v1-vSc_initial);
deltav2 = norm(v2-vSc_final);
deltavTot = deltav1 + deltav2;
