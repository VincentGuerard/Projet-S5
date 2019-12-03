%Compensateur calculer (Angles et hauteur z)
numCompAngles = 1.0e09*[0.000230814504577 0.026170635469275 0.897433010674818 9.699372248008062];
denCompAngles = 1.0e07*[0.000000100000000 0.000823402569667 1.694979479335234 0];
numCompz = 1.0e11*[0.000031609546854 0.004304100455063 0.142047621936623 1.354942422178642];
denCompz = 1.0e05*[0.000010000000000 0.017380342819025 5.055406372627656 0];
CompAngles = tf(numCompAngles, denCompAngles);
Compz = tf(numCompz, denCompz);

%Continue à Discret (s à z)
testdiscret(CompAngles)
testdiscret(Compz)

%Le lieu de bode est semblable pour les compensateurs!!