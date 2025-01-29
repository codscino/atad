%% function to get rf with dt instead of dtheta
function [rf, vf] = FGKepler_dt2(r0,v0,dt,mu)
    
    R0 = norm(r0);
    V0 = norm(v0);
    a = mu / ((2*mu)/R0 - V0^2);

    % safe check
    if a < 0
        a
        warning('DOH! a is hyperbolic, the method does not converge!')
        %v0 = [NaN NaN NaN];
        rf = [NaN NaN NaN];
        vf = [NaN NaN NaN];
        return
    end
    
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
    
    % II part
    Rf = norm(rf);
    dG = 1 - (a/Rf)*(1-cos(dE));
    vf = 1/G * (dG*rf - r0);

end
