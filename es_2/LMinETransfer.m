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