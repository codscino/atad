%% initial data
hp = 400; %km
ha = 6000;
dt = 0.65*3600;

R = getAstroConstants('Earth', 'Radius');
rp = hp + R;
ra = ha + R;

%% function

function [E,niteration] = newton_raphson(rp,ra,mu,dt,tol)

    a = 0.5*(rp+ra);
    e = (ra-rp)/(ra+rp);
    n = sqrt(mu/a^3);
    E0 = 0; %periapsis
    M0 = 0; %periapsis
    M1 = n*dt;

    err = 10;
    f1 = (M1-e*sin(M1))-M1;
    df1 = 1-e*cos(M1);
    Ek = M1 - f1/df1;
    niteration = 0;

    while err>tol
        niteration = niteration +1;
        fEk = (Ek-e*sin(Ek))- M1;
        dfEk = 1-e*cos(Ek);
        err = abs(fEk);
        Ek = Ek - (fEk/dfEk);    
    end
   
    E = Ek;
end


%% final results
mu = getAstroConstants('Earth', 'Mu');
tol = 10^(-6);
[E,niteration] = newton_raphson(rp,ra,mu,dt,tol)