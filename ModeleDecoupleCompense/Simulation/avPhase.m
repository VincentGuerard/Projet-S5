function compensateur = avPhase(poles_desire, FTBO, ajoutPhase, compensateur_simple)
%compensateur = avPhase(poles_desire, FTBO, ajoutPhase, compensateur_simple)
%   Detailed explanation goes here

[num den] = tfdata(FTBO, 'v');
Gspd = polyval(num, poles_desire(1))/polyval(den, poles_desire(1));
phaseGs = -(2*pi - angle(Gspd));
deltaPhase = -pi - phaseGs + ajoutPhase;

if deltaPhase < (65/360)*2*pi || compensateur_simple
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)/(poles_desire(1)-p))*abs(Gspd));

    Ga = tf(Ka*[1 -z], [1 -p]);
else
    deltaPhase = deltaPhase / 2;
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)^2/(poles_desire(1)-p)^2)*abs(Gspd));
    Ka = sqrt(Ka);

    Ga = tf(Ka*[1 -z], [1 -p]);
    Ga = series(Ga, Ga);
end

compensateur = Ga;

end