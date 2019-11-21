function compensateur = avPhase(poles_desire, FTBO, ajoutPhase, compensateur_simple)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[num den] = tfdata(FTBO, 'v');
Gspd = polyval(num, poles_desire(1))/polyval(den, poles_desire(1));
phaseGs = -(2*pi - angle(Gspd));
deltaPhase = -pi - phaseGs + ajoutPhase

if deltaPhase < (60/360)*2*pi || compensateur_simple
    % Compensateur simple
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)/(poles_desire(1)-p))*abs(Gspd));

    Ga = tf(Ka*[1 -z], [1 -p]);
    
elseif deltaPhase > (60/360)*2*pi && deltaPhase < (120/360)*2*pi
    % Compensateur double
    deltaPhase = deltaPhase / 2;
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)/(poles_desire(1)-p))*abs(Gspd))
    Ka = sqrt(Ka);

    Ga = tf(Ka*[1 -z], [1 -p]);
    Ga = Ga * Ga;

elseif deltaPhase > (120/360)*2*pi && deltaPhase < (180/360)*2*pi
    % Compensateur triple
    deltaPhase = deltaPhase / 3;
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)/(poles_desire(1)-p))*abs(Gspd))
    Ka = sqrt(Ka);

    Ga = tf(Ka*[1 -z], [1 -p]);
    Ga = Ga * Ga * Ga;
    
elseif deltaPhase > (180/360)*2*pi
    % Compensateur quadruple
    deltaPhase = deltaPhase / 4;
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)/(poles_desire(1)-p))*abs(Gspd))
    Ka = sqrt(Ka);

    Ga = tf(Ka*[1 -z], [1 -p]);
    Ga = Ga * Ga * Ga * Ga;
end

compensateur = Ga;

end