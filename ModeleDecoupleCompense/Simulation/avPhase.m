function compensateur = avPhase(poles_desire, FTBO, ajoutPhase, force_compensateur)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[num den] = tfdata(FTBO, 'v');
Gspd = polyval(num, poles_desire(1))/polyval(den, poles_desire(1));
phaseGs = -(2*pi - angle(Gspd));
deltaPhase = -pi - phaseGs + ajoutPhase;
rad2deg(deltaPhase)

if deltaPhase < (60/360)*2*pi || force_compensateur == 1
    % Compensateur simple
    disp('Compensateur simple')
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)/(poles_desire(1)-p))*abs(Gspd));

    Ga = tf(Ka*[1 -z], [1 -p]);
    
elseif deltaPhase > (60/360)*2*pi && deltaPhase < (120/360)*2*pi || force_compensateur == 2
    % Compensateur double
    disp('Compensateur double')
    deltaPhase = deltaPhase / 2;
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)^2/(poles_desire(1)-p)^2)*abs(Gspd))
    Ka = sqrt(Ka);

    Ga = tf(Ka*[1 -z], [1 -p]);
    Ga = Ga * Ga;

elseif deltaPhase > (120/360)*2*pi && deltaPhase < (180/360)*2*pi || force_compensateur == 3
    % Compensateur triple
    disp('Compensateur triple')
    deltaPhase = deltaPhase / 3;
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)^3/(poles_desire(1)-p)^3)*abs(Gspd))
    Ka = nthroot(Ka,3);

    Ga = tf(Ka*[1 -z], [1 -p]);
    Ga = Ga * Ga * Ga;
    
elseif deltaPhase > (180/360)*2*pi && deltaPhase < (240/360)*2*pi  || force_compensateur == 4
    % Compensateur quadruple
    disp('Compensateur quadruple')
    deltaPhase = deltaPhase / 4;
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)^4/(poles_desire(1)-p)^4)*abs(Gspd))
    Ka = nthroot(Ka,4);

    Ga = tf(Ka*[1 -z], [1 -p]);
    Ga = Ga * Ga * Ga * Ga;

elseif deltaPhase > (240/360)*2*pi && deltaPhase < (300/360)*2*pi || force_compensateur == 5
    % Compensateur quadruple
    disp('Compensateur cinqtuple')
    deltaPhase = deltaPhase / 5;
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)^5/(poles_desire(1)-p)^5)*abs(Gspd))
    Ka = nthroot(Ka,5);

    Ga = tf(Ka*[1 -z], [1 -p]);
    Ga = Ga * Ga * Ga * Ga * Ga;
    
elseif deltaPhase > (300/360)*2*pi || force_compensateur == 6
    % Compensateur sixtuple
    disp('Compensateur sixtuple')
    deltaPhase = deltaPhase / 6;
    alpha = pi - (pi - angle(poles_desire(1)));
    phiz = (alpha + deltaPhase)/2;
    phip = (alpha - deltaPhase)/2;

    p = real(poles_desire(1)) - imag(poles_desire(1))/tan(phip);
    z = real(poles_desire(1)) - imag(poles_desire(1))/tan(phiz);

    % Gain Ka
    Ka = 1/abs(((poles_desire(1)-z)^6/(poles_desire(1)-p)^6)*abs(Gspd))
    Ka = nthroot(Ka,6);

    Ga = tf(Ka*[1 -z], [1 -p]);
    Ga = Ga * Ga * Ga * Ga * Ga * Ga;
end

compensateur = Ga;

end
